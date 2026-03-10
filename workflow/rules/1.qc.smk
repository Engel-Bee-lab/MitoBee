"""
Rules for quality control and quality assurance - Illumina paired end reads 
"""
import glob
import os

#quality control rules here
rule fastp:
    input:
        r1 = lambda wc: os.path.join(input_dir, f"{wc.sample}{pattern_r1}.{extn}"),
        r2 = lambda wc: os.path.join(input_dir, f"{wc.sample}{pattern_r2}.{extn}")
    output:
        r1 = os.path.join(dir_fastp,"{sample}_R1.fastq.gz"),
        r2 = os.path.join(dir_fastp,"{sample}_R2.fastq.gz"),
        stats = os.path.join(dir_fastp,"{sample}.stats.json"),
        html = os.path.join(dir_fastp,"{sample}.stats.html")
    conda:
        os.path.join(dir_env, "fastp.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        # Default params in this command 
        #   --qualified_quality_phred = 15 
        #   --unqualified_percent_limit 40 
        #   --n_base_limit 5 
        #   --length_required 30
        
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
            --thread {threads} \
            --json {output.stats} \
            --html {output.html}
        """
    
