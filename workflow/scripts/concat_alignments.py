#!/usr/bin/env python

import os
from Bio import SeqIO

# ---------------------------
# Snakemake variables
# ---------------------------
# snakemake.input.tmp_dir -> folder with *_aligned.faa files
# snakemake.output.fasta -> concatenated FASTA
# snakemake.output.partitions_txt -> Partitions.txt
# snakemake.output.partitions_nex -> Partitions.nex

try:
    tmp_dir = os.path.abspath(snakemake.params.tmp_dir)
    output_fasta = os.path.abspath(snakemake.output.fasta)
    partition_txt = os.path.abspath(snakemake.output.partitions_txt)
    partition_nex = os.path.abspath(snakemake.output.partitions_nex)
except NameError:
    # Mock values for testing outside Snakemake
    tmp_dir = os.path.abspath("input_tmp_dir")
    output_fasta = os.path.abspath("output_concat.fasta")
    partition_txt = os.path.abspath("output_partitions.txt")
    partition_nex = os.path.abspath("output_partitions.nex")

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
    if not os.path.exists(aln_file):
        print(f"[INFO] {gene}_aligned.faa not found, skipping")
        continue

    records = list(SeqIO.parse(aln_file, "fasta"))
    gene_lengths[gene] = len(records[0].seq)

    for rec in records:
        gid = rec.id
        if gid not in genomes:
            genomes[gid] = []
        genomes[gid].append(str(rec.seq))

# ---------------------------
# Write concatenated FASTA
# ---------------------------
with open(output_fasta, "w") as outfh:
    for gid, seqs in genomes.items():
        concatenated = ''.join(seqs)
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