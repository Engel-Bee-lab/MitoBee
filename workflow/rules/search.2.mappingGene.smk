rule host_mapping_search_gene:
    input:
        r1 = os.path.join(dir_fastp,"{sample}_R1.fastq.gz"),
        r2 = os.path.join(dir_fastp,"{sample}_R2.fastq.gz"),
        host=config['args']['ref_set']
    output:
        all_bam=os.path.join(dir_hostsearch,"{sample}_temp_gene.bam"),
    params:
        host_group= os.path.join(dir_hostsearch, "ref_genes.fasta"),
        dirs=os.path.join(dir_hostsearch),
    conda:
        os.path.join(dir_env, "bowtie2.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        set -euo pipefail
        rm -rf {params.host_group}
        cat {input.host}/*.fa* >> {params.host_group}
        
        #moving to bowtie2 since its developed for shorter references too 
        #if I use minimap2 then I need to decrease the kmer and -w params, based on the gene length
        #bowtie2, I dont need to change this
        
        #index reference genes 
        bowtie2-build {params.host_group} {params.dirs}/gene_index
        bowtie2 -x {params.dirs}/gene_index -1 {input.r1} -2 {input.r2} -S {output.all_bam} | \
            samtools sort -@ {threads} -o {output.all_bam} -

        samtools index {output.all_bam}
        """

"""
Given the mapping metrics considered should be based on if the reference set are 
- within the same species 
- within the same genus
- across genera
"""
rule host_mapping_metrics_gene:
    input:
        bam=os.path.join(dir_hostsearch,"{sample}_temp_gene.bam")
    output:
        metrics=os.path.join(dir_hostsearch, "{sample}_all_idxstats_gene.txt"),
        primary=os.path.join(dir_hostsearch,"{sample}_primary_idxstats_gene.txt"),
        coverage=os.path.join(dir_hostsearch,"{sample}_primary_coverage_gene.txt"),
        strict=os.path.join(dir_hostsearch,"{sample}_strict_idxstats_gene.txt")
    params:
        primary_bam=os.path.join(dir_hostsearch,"{sample}_primary_gene.bam"),
        strict_bam=os.path.join(dir_hostsearch,"{sample}_strict_gene.bam")
    conda:
        os.path.join(dir_env, "minimap2.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        set -euo pipefail

        #linient mapping metrics:
        #outside genus
        samtools idxstats {input.bam} > {output.metrics}
        
        #primary alignments only (remove secondary and supplementary alignments):
        #within genus
        samtools view  -b -F 0x900 {input.bam} > {params.primary_bam}
        samtools sort -@ {threads} -o {params.primary_bam}.sorted.bam {params.primary_bam}
        mv {params.primary_bam}.sorted.bam {params.primary_bam}
        samtools index {params.primary_bam}
        samtools idxstats {params.primary_bam} > {output.primary}
        samtools coverage {params.primary_bam} > {output.coverage}

        #strict mappimng metrics (conservative metrics with only primary alignments and high mapping quality):
        #within subspecies
        samtools view  -b -F 0x900 -q 30 {input.bam} > {params.strict_bam}
        samtools sort -@ {threads} -o {params.strict_bam}.sorted.bam {params.strict_bam}
        mv {params.strict_bam}.sorted.bam {params.strict_bam}
        samtools index {params.strict_bam}
        samtools idxstats {params.strict_bam} > {output.strict}

        rm {params.primary_bam} {params.strict_bam}
        """

"""
Now adding a rule for scoring mechanism to determine the best mapping reference for each sample.
"""
rule host_mapping_score_gene:
    input:
        primary_idx = os.path.join(dir_hostsearch,"{sample}_primary_idxstats_gene.txt"),
        primary_cov = os.path.join(dir_hostsearch,"{sample}_primary_coverage_gene.txt"),
        strict_idx = os.path.join(dir_hostsearch,"{sample}_strict_idxstats_gene.txt")
    output:
        summary = os.path.join(dir_reports, "Reference_search", "{sample}_host_ranking_gene.tsv")
    localrule: True
    script:
        "../scripts/score_hosts.py"
