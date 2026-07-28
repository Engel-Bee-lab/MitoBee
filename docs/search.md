# Search to determine the most closely related reference gene/genome to use  

This module helps identify an appropriate mitochondrial reference genome to use with the `mitobee run` workflow. It is particularly useful when a suitable reference mitogenome is not already known for the organism of interest.

The workflow evaluates how well sequencing reads map to a set of candidate mitochondrial references and ranks them based on mapping statistics. The highest-ranked reference genome can then be used for downstream mitogenome reconstruction.

### Method used  
Reads are mapped against the provided reference set, and several mapping metrics are calculated for each reference sequence:
- Total number of mapped reads
- Number of high-quality (unique) alignments
- Breadth of coverage
- Mean read depth

## Inputs
Note: There is likely no reason to run all the metagenomes here, maybe the subset of the samples, based on the host diversity. 

- Metagenomic reads, paired end saved in one directory. So the folder structure should be as following

        ├── Metagenomes
        │   ├── meta1_R1.fastq
        │   ├── meta1_R2.fastq
        │   ├── meta2_R1.fastq
        │   ├── meta2_R2.fastq
        │   ├── meta3_R1.fastq
        │   ├── meta3_R2.fastq
        │   ├── ...

    The file extension can be .fastq, .fastq.gz, .fq etc. It doesnt matter as long as the correct ending is added to the command 
    Also, the forward and the reverse reads can be '_R1'/'_R2' or '_1'/'_2, as long as they are defined in the command.

- Reference datasets: Includes either mitogenomes or specific genes
  Based on this, add params `--mode mitogenome` or `--mode genes`.

## Example command 

```
    #If a closely related mitochondrial genome is not available, but a gene is, like cox  or rRNA genes
    #Download the reference genes you would like to use of the closely related genomes

    #to search against mitogenomes refernece set
    mitobee gene --input test-files/mitogenomes --extn fasta --ref_seq  test-files/ref --output output -k all --mode mitogenome

    #to search against mitogenomes refernece gene set 
    mitobee gene --input test-files/mitogenomes --extn fasta --ref_seq  test-files/ref --output output -k all --mode genes 
```

## Outputs 
For each sample, the workflow produces intermediate mapping statistics as well as a final ranked summary:
- PROCESSING 
- REPORTS