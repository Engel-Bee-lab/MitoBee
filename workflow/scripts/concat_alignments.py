#!/usr/bin/env python

import os
from Bio import SeqIO
import argparse

def get_paths():
    # Standalone execution context
    parser = argparse.ArgumentParser()
    parser.add_argument("--tmp-dir", required=True)
    parser.add_argument("--output-fasta", required=True)
    parser.add_argument("--partition-txt", required=True)
    parser.add_argument("--partition-nex", required=True)
    args = parser.parse_args()

    return (
        os.path.abspath(args.tmp_dir),
        os.path.abspath(args.output_fasta),
        os.path.abspath(args.partition_txt),
        os.path.abspath(args.partition_nex),
    )

tmp_dir, output_fasta, partition_txt, partition_nex = get_paths()

os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
os.makedirs(os.path.dirname(partition_txt), exist_ok=True)

# ---------------------------
# Fixed mitochondrial protein gene order
# ---------------------------
genes = [
    "atp6", "atp8", "cob", "cox1", "cox2", "cox3",
    "nad1", "nad2", "nad3", "nad4", "nad4L", "nad5", "nad6"
]

# ---------------------------
# Collect sequences per genome
# ---------------------------
genomes = {}         # genome_id -> list of sequences in order
gene_lengths = {}    # gene -> length of aligned sequence

for gene in genes:
    aln_file = os.path.join(tmp_dir, f"{gene}_aligned.faa")
    print(aln_file)
    if not os.path.exists(aln_file):
        print(f"[INFO] {gene}_aligned.faa not found, skipping")
        continue

    records = list(SeqIO.parse(aln_file, "fasta"))
    if not records:
        print(f"[WARN] {gene}_aligned.faa is empty, skipping")
        continue
    
    # Record gene length from first sequence
    gene_lengths[gene] = len(records[0].seq)

    for rec in records:
        gid = rec.id
        seq = str(rec.seq)
        if gid not in genomes:
            genomes[gid] = seq
        else:
            genomes[gid] += seq

# ---------------------------
# Write concatenated FASTA
# ---------------------------
with open(output_fasta, "w") as outfh:
    for gid, concatenated in filtered_genomes.items():
        outfh.write(f">{gid}\n{concatenated}\n")

print(f"[INFO] Concatenated alignment written to {output_fasta}")

# ---------------------------
# Write Partitions.txt
# ---------------------------
curr_start = 1
with open(partition_txt, "w") as f:
    for gene in genes:
        if gene not in gene_lengths:
            continue
        length = gene_lengths[gene]
        curr_stop = curr_start + length - 1
        f.write(f"AA, {gene} = {curr_start}-{curr_stop}\n")
        curr_start = curr_stop + 1

# ---------------------------
# Write Partitions.nex
# ---------------------------
with open(partition_txt) as f:
    lines = [line.strip() for line in f if line.strip()]

with open(partition_nex, "w") as f:
    f.write("#NEXUS\nbegin sets;\n")
    for line in lines:
        parts = line.split(',')
        if len(parts) != 2:
            continue
        _, rest = parts
        name, coords = rest.split('=')
        name = name.strip().replace(" ", "_")
        coords = coords.strip()
        f.write(f"  charset {name} = {coords};\n")
    f.write("end;\n")

print(f"[INFO] Partitions written to {partition_txt} and {partition_nex}")