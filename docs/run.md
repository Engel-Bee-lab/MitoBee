# Reference-based Mitogenome reconstruction in metagenome

This module of **MitoBee** reconstructs a mitochondrial genome by extracting reads that map to a user-provided reference mitogenome. The goal is to recover a consensus sequence for each sample using a closely related mitochondrial reference. Please use the search option 

## Method overview

1. **Quality control**
   FastP for quality control

2. **Reference Mapping**
   Trimmed reads are mapped against a reference mitochondrial genome supplied by the user. The reference should ideally come from a closely related species to maximize mapping accuracy.

3. **Read Extraction and Consensus Generation**
   Reads that successfully map to the reference are used to build a consensus mitochondrial genome for each sample.

4. **Low-Coverage Masking**
   Positions with insufficient read support are masked to prevent unreliable base calls:

   * Sites with **read depth ≤ 10** are replaced with **N** in the consensus sequence.

5. **Quality Filtering of Consensus Genomes**
   To ensure only reliable assemblies are reported, each consensus genome is evaluated based on the proportion of masked positions.

   * If **more than 33% of the consensus sequence consists of Ns**, the assembly is considered **low quality**.
   * In this case, the mitogenome for that sample **fails quality control and is not reported**.

## Input 
- Metagenomic reads, paired end saved in one directory. So the folder structure should be as following
      
      ├── Metagenomes
      │   ├── meta1_R1.fastq
      │   ├── meta1_R2.fastq
      │   ├── meta2_R1.fastq
      │   ├── meta2_R2.fastq
      │   ├── meta3_R1.fastq
      │   ├── meta3_R2.fastq
      │   ├── ....
  
  The file extension can be .fastq, .fastq.gz, .fq etc. It doesnt matter as long as the correct ending is added to the command 
  Also, the forward and the reverse reads can be '_R1'/'_R2' or '_1'/'_2, as long as they are defined in the command. 

- Reference host
  Provide the fasta file containing mitochondrial genome, preferably complete. 


## Unsure of which reference gene/genome to use
If you have sampled from a unique reference host, whose chromsome or mitochondrial genome still hasn't been sequenced. Then running below steps presents a challenge.
Since the below `run` workflow is best supported when used with a closely related genome. To determine which available reference genome is the best suited to use, run 
[`mitobee search`](https://mitobee.readthedocs.io/en/latest/search/) command usage. 

## Example command 
There are two modes for running this command
1. **mitogenome mode**: Where the complete mtDNA will be used to build a consensus mtDNA from each metagenome sample

```
    #Running mitobee with test files available in the repo
    mitobee run --input <path to metagenome directory> --extn fastq.gz \
         --pattern_r1 _R1 --pattern_r2 _R2 \
         --host_seq <host reference mtgenome> \
         --output output --mode mitogenome
```

2. **gene mode**: In this case, the user can provide a bunch of clsoely
Reasoning: Often the CO1 genes or the cytB/cob gene are sequenced for more hosts, but not the entire mitogenome. In these cases, the user can also run this workflow to build consensus genes from metagenomes. 

```
    #Running mitobee with test files available in the repo
    mitobee run --input <path to metagenome directory>  --extn fastq.gz \
         --pattern_r1 _R1 --pattern_r2 _R2 \
         --host_seq <reference host gene> \
         --output output --mode gene
```

## Output

For each sample, the module produces:
- PROCESSING: All the intermediate results from all the tools run
  - 1_fastp: FastP outputs for quality control of metagenomes 
  - Host_cleaned: mapping results, bam files, fastq files including reads mapped
  - Host_cleaned/mitogenome: Building the consensus mitogenome intermediate files

- REPORTS:
  - QC/<sample>_qc_report.txt: Includes the QC stats for each metagenome. For example:  

      Sample: test1 \
      Total Reads: 12768112 \ 
      QC Reads: 12768114 \
      Mapped Reads: 3447 \
      Percentage of Mapped Reads: 0.03% \

  - final_host_qc_summary.txt: All the QC stats in one file 

    |sample|	total_reads|	QC_reads|	mapped_reads|	percent_mapped|
    |------|------|-----|-----|-----|
    |tSTS|	12768112|	12768114|	3447|	0.03|
    |cdbjs|	362050|	362050|	498	|0.14|
    |test1|	6655074|	6659740|	3631293|	54.53|
    |test2|	5394142|	5394269|	376121|	6.97|

  - mitogenome_consensus_summary.tsv: Consensus mitogenome stats

      |filename|	header|	length|	GC_content|	N_count|	N_fraction|	QC_status|
      |test1_consensus.fasta|	test1_CM009947.2|	16471|	15.03|	139|	0.0084|	PASS|
      |test2_consensus.fasta|	test2_CM009947.2|	16471| 14.95|	216|	0.0131|	PASS|

  - mitogenome/<sample name>_consensus.fasta: Includes the resulting mitogenome generated from the mapped reads from the metagenome. If they pass the quality. 

  Samples that fail the masking threshold (>33% Ns) will **not produce a final mitogenome sequence**, indicating that insufficient mitochondrial signal was present in the dataset.

## Recommendations

* Use a **closely related mitochondrial reference genome** to improve mapping efficiency and consensus accuracy.
* Samples with **very low mitochondrial read content** may fail the quality threshold due to extensive masking.
