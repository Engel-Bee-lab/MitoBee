#!/bin/bash

#SBATCH --account pengel_default
#SBATCH --job-name mitomee_dev
#SBATCH --output dev-%j.error
#SBATCH --error dev-%j.out
#SBATCH -N 1
#SBATCH --cpus-per-task 16
#SBATCH --partition=cpu
#SBATCH --mem=150G
#SBATCH --time=8:00:00


#mitobee run --input test-files/metagenomes --extn fastq.gz \
#    --pattern_r1 _R1 --pattern_r2 _R2 --cores 16 \
#     --host_seq test-files/am-dh4.fasta --profile slurm \
#     --conda-frontend mamba --output output 

#build tree with host seq, assembled mitogenomes. 
#Also add other reference mitogenomes to  output/REPORTS/mitogenomes. No need to add the reference sequence, it will be included
#Then run the subcommand to build tree 

mitobee tree --input test-files/mitogenomes --extn fasta --output output -k all --profile slurm --conda-frontend mamba --cores 16

#mitobee search --input test-files/metagenomes/ --extn fastq.gz --pattern_r1 _R1 --pattern_r2 _R2 --ref_set test-files/ref-set-genome/ --output output --mode mitogenome -k --profile slurm --conda-frontend mamba