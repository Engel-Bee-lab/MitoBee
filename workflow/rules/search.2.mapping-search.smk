"""
Rules for host searchand mapping. 
The goal of this rule is to search for the host sequence that best maps to the ref set
"""

rule host_mapping_search:
    input:
        r1 = os.path.join(dir_fastp,"{sample}_R1.fastq.gz"),
        r2 = os.path.join(dir_fastp,"{sample}_R2.fastq.gz"),
        host=config['args']['ref_set']
    output:
        all_bam=os.path.join(dir_hostsearch,"{sample}_temp.bam"),
    params:
        host_group= os.path.join(dir_hostsearch, "ref_group.fasta"),
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
        rm -rf {params.host_group}
        cat {input.host}/*.fasta >> {params.host_group}

        minimap2 -ax sr -t {threads} {params.host_group} {input.r1} {input.r2} \
            | samtools sort -@ {threads} -o {output.all_bam} -

        samtools index {output.all_bam}
        """

"""
Given the mapping metrics considered should be based on if the reference set are 
- within the same species 
- within the same genus
- across genera
"""
rule host_mapping_metrics:
    input:
        bam=os.path.join(dir_hostsearch,"{sample}_temp.bam")
    output:
        metrics=os.path.join(dir_hostsearch, "{sample}_all_idxstats.txt"),
        primary=os.path.join(dir_hostsearch,"{sample}_primary_idxstats.txt"),
        coverage=os.path.join(dir_hostsearch,"{sample}_primary_coverage.txt"),
        strict=os.path.join(dir_hostsearch,"{sample}_strict_idxstats.txt")
    params:
        primary_bam=os.path.join(dir_hostsearch,"{sample}_primary.bam"),
        strict_bam=os.path.join(dir_hostsearch,"{sample}_strict.bam")
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
rule host_mapping_score:
    input:
        primary_idx = os.path.join(dir_hostsearch,"{sample}_primary_idxstats.txt"),
        primary_cov = os.path.join(dir_hostsearch,"{sample}_primary_coverage.txt"),
        strict_idx = os.path.join(dir_hostsearch,"{sample}_strict_idxstats.txt")
    output:
        summary = os.path.join(dir_hostsearch,"{sample}_host_ranking.tsv")
    localrule: True
    script:
        "../scripts/score_hosts.py"