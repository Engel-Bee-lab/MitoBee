"""
Rules for quality control and quality assurance - Illumina paired end reads 
"""

#quality control rules here
rule fastp:
    input:
        r1 = os.path.join(input_dir, PATTERN_R1),
        r2 = os.path.join(input_dir, PATTERN_R2)
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
        
        if [ -f {output.r1} ] && [ -f {output.r2} ]; then
            echo "fastp output files already exist. Skipping..."
        else
            fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
                --thread {threads} \
                --json {output.stats} \
                --html {output.html}
        fi
        
        if [ ! -s {output.r1} ] || [ ! -s {output.r2} ]; then
            echo "Error: One of the output files is empty. Please check the input files and parameters."
            touch {output.r1} {output.r2} # create empty files to avoid snakemake errors
        fi
        """
    
