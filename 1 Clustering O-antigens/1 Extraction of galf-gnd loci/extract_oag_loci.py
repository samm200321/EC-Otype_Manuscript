import pandas as pd
import matplotlib.pyplot as plt
import os
from Bio import SeqIO
import gzip
import subprocess
from tqdm import tqdm
from Bio.SeqRecord import SeqRecord
from pathlib import Path

#to run download assemblies to folder called ec_refseq
fasta_dir = Path("ec_refseq")

df = [f for f in os.listdir('ec_refseq') if f.endswith('.fna.gz')]
#df.set_index("Name", inplace=True)

print(df)



#searchfor = ["wzy", "wzx"]
#df_filtered = df[df["GeneScores"].str.contains('|'.join(searchfor))]
df_filtered = df
#df_filtered.to_csv('filtered_genomes.tsv', sep='\t', index=False)


#counts = df_filtered['O-type'].value_counts()

gnd = SeqIO.read("gnd.fasta", "fasta")
galF = SeqIO.read("galF.fasta", "fasta")


print(gnd)
print(galF)

#print(df_filtered.iloc[0].name)

#record = SeqIO.read(f"fna/{df_filtered.iloc[10].name}.fna", "fasta")




def blast(combination, pid, ev):
    #cmd = f"blastn -query \"{combination[0]}\" -subject {combination[1]} -out \"tmp/{combination[1]}.csv\""
    cmd = f"blastn -query \"{combination[0]}\" -subject <(gzip -dc {combination[1]}) -out \"tmp/{combination[1]}.csv\""
    # print(cmd)
    cmd += f" -outfmt 6"
    subprocess.run(["bash", "-c", cmd])

    if os.path.getsize(f"tmp/{combination[1]}.csv") == 0:
        return None

    return pd.read_csv(f"tmp/{combination[1]}.csv", delimiter="\t", header=None)

o_loc_db = []
o_loc_db_long = []
not_same_contig = 0
no_border = 0
for i in tqdm(range(0, len(df_filtered))):


    galF_start = blast(["galF.fasta", f"ec_refseq/{df_filtered[i]}"], 60, 0.01)

    gnd_start = blast(["gnd.fasta", f"ec_refseq/{df_filtered[i]}"], 60, 0.01)
    if galF_start is None or gnd_start is None:
        print("galf/gnd not found\n")
        no_border += 1
        #print(df_filtered.iloc[i]["O-type"])
        continue
    if galF_start.iloc[0, 1] != gnd_start.iloc[0,1]:
        not_same_contig += 1
        continue

    record = None
    target_contig = galF_start.iloc[0, 1]
    with gzip.open(f"ec_refseq/{df_filtered[i]}", "rt") as handle:
        print(handle)
        for tmp_record in SeqIO.parse(handle, "fasta"):
            if tmp_record.id == target_contig:
                record = tmp_record
                continue


    o_loc_seq = ""
    if len(galF_start) > 1 or len(gnd_start) > 1:
        print("more than 1 hit")
    if galF_start.iloc[0, 8]-gnd_start.iloc[0, 9] > 0:

        o_loc_seq = record[gnd_start.iloc[0, 8]:galF_start.iloc[0, 9]].seq
    else:

        o_loc_seq = record[galF_start.iloc[0, 9]:gnd_start.iloc[0, 8]].seq.reverse_complement()

    o_loc = SeqRecord(
        o_loc_seq,
        id=record.id,
        name=f"O-antigen  Locus",
    )
    if len(o_loc.seq) > 30000:
        o_loc_db_long.append(o_loc)
    else:


        o_loc_db.append(o_loc)


SeqIO.write(o_loc_db, "O-antigen_db.fasta_all.fasta", "fasta")
SeqIO.write(o_loc_db_long, "O-antigen_db.fasta_long_all.fasta", "fasta")
print(f"{not_same_contig} not on the same contig")
print(f"{no_border} no galF/gnd")



