"""
Extra steps with the host reads 
    - map to the the reference mitogenome
"""
from glob import glob

rule bam_sort:
    input:
        mapped_bam = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mapped.bam"),
    output:
        sort_bam=os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mapped.sorted.bam"),
        depth=os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mitogenome_snps_depth.txt")
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
        if [ -f {output.sort_bam} ]; then
            echo "Output file found, so this run looks liks its run. Skipping..."
            exit 0
        else
        
            #prep the bam
            samtools sort -o {output.sort_bam} {input.mapped_bam}
            samtools index {output.sort_bam}
            samtools depth {output.sort_bam} > {output.depth}
        fi
        """

rule host_mito_snps:
    input:
        sort_bam = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mapped.sorted.bam"),
        host= config['args']['host_seq']
    output:
        vcf = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mitogenome_snps.vcf.gz"),
        filterred_vcf = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mitogenome_snps.filtered.vcf.gz")
    params:
        prefix = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mitogenome"),
        stats=os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mitogenome_snps.stats.txt"),
    conda:
        os.path.join(dir_env, "bcftools.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        set -euo pipefail
        if [ -f {output.vcf} ]; then
            echo "Output vcf file found, so this run looks liks its run. Skipping..."
            exit 0
        else

            #call variants
            bcftools mpileup -f {input.host} -Q 20 -q 20 {input.sort_bam} | \
                bcftools call --ploidy 1 -mv -Ov -o {output.vcf}
            bcftools index {output.vcf}

            #filter the SNPs, QUAL>30 means 0.1% error rate, DP>10 means at least 10 reads support
            bcftools filter -i 'QUAL>30 && DP>10' {output.vcf} -Oz -o {output.filterred_vcf}
            bcftools index {output.filterred_vcf}

            #getting the stats
            bcftools stats {output.filterred_vcf} > {params.stats}
        fi
        """

rule normalise_vcfs:
    input:
        filterred_vcf = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mitogenome_snps.filtered.vcf.gz")
    output:
        norm_vcf = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mitogenome_snps.filtered.norm.vcf.gz")
    params:
        host= config['args']['host_seq'],
        temp=os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mitogenome_snps.filtered.temp.vcf.gz"), 
    conda:
        os.path.join(dir_env, "bcftools.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        set -euo pipefail
        if [ -f {output.norm_vcf} ]; then
            echo "Normalized VCF already exists. Skipping..."
            exit 0
        else
            bcftools norm -f {params.host} -m -any {input.filterred_vcf} -Oz -o {params.temp} 
            bcftools index {params.temp}

            #then keep only SNPs
            bcftools view -v snps {params.temp} -Oz -o {output.norm_vcf}
            bcftools index {output.norm_vcf}
        fi
        """

rule generate_allele_frequency:
    input:
        norm_vcf = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mitogenome_snps.filtered.norm.vcf.gz")
    output:
        af_table = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mitogenome_snps.allele_frequency.txt")
    params:
        prefix = "{sample}",
    conda:
        os.path.join(dir_env, "bcftools.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        set -euo pipefail
        if [ -f {output.af_table} ]; then
            echo "Allele frequency table already exists. Skipping..."
            exit 0
        else
            bcftools query -f '%POS\t[%AD]\n' "{input.norm_vcf}" |   awk -v s={params.prefix} '{{ 
                split($2,a,",");
                if (length(a) == 2) {{
                    total = a[1] + a[2];
                    if (total > 0)
                        print s "\t" $1 "\t" a[2]/total
                }}
            }}' > {output.af_table}
        fi
        """

rule merge_vcf:
    input:
        expand(os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mitogenome_snps.filtered.norm.vcf.gz"), sample=sample_names)
    output:
        merged_vcf = os.path.join(dir_hostcleaned, "mitogenome", "merged_mitogenome_snps.filtered.norm.vcf.gz")
    params:
        temp=os.path.join(dir_hostcleaned, "mitogenome", "merged_mitogenome_snps.temp.vcf.gz"),
        folder=os.path.join(dir_hostcleaned, "mitogenome")
    conda:
        os.path.join(dir_env, "bcftools.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        set -euo pipefail
        if [ -f {output.merged_vcf} ]; then
            echo "Merged VCF already exists. Skipping..."
            exit 0
        else
            bcftools merge -Oz -o {params.temp} {params.folder}/*_mitogenome_snps.filtered.norm.vcf.gz

            mv {params.temp} {output.merged_vcf}
            bcftools index {output.merged_vcf}
        fi
        """

rule low_coverage_bed:
    input:
        sort_bam = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_mapped.sorted.bam")
    output:
        low_cov_bed = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_low_cov.bed")
    conda:
        os.path.join(dir_env, "minimap2.yaml") 
    params:
        min_depth = 10 #mask regions with less than 10x coverage
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        set -euo pipefail
        # -a outputs all positions (including zero depth)
        # We filter for positions where depth ($3) is LESS than the threshold
        samtools depth -a {input.sort_bam} \
            | awk -v min={params.min_depth} '$3 < min {{print $1"\t"$2-1"\t"$2}}' > {output.low_cov_bed}
        """

rule snp_alignment:
    input:
        merged_vcf = os.path.join(dir_hostcleaned, "mitogenome", "merged_mitogenome_snps.filtered.norm.vcf.gz"),
        host= config['args']['host_seq'],
        low_cov_bed = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_low_cov.bed")
    output:
        consensus_fasta = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_consensus.fasta")
    params:
        sample = "{sample}",
    conda:
        os.path.join(dir_env, "bcftools.yaml")
    params:
        sample = "{sample}",
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        set -euo pipefail
        if [ -f {output.consensus_fasta} ]; then
            echo "SNP alignment fasta already exists. Skipping..."
            exit 0
        else
            SAMPLE_FULL=$(bcftools query -l {input.merged_vcf} | grep "{params.sample}")

            # Generate consensus masking zero coverage
            bcftools consensus -s "$SAMPLE_FULL" -f {input.host} -m {input.low_cov_bed} {input.merged_vcf} > {output.consensus_fasta}

            # Add sample name to fasta header
            sample="{params.sample}"
            sed -i "s/>/>${{sample}}_/g" {output.consensus_fasta}
        fi
        """

rule qc_consensus:
    input:
        fasta = os.path.join(dir_hostcleaned, "mitogenome", "{sample}_consensus.fasta")
    output:
        output_fasta = os.path.join(dir_reports, "mitogenome", "{sample}_done.txt"),
    params:
        max_frac = 0.333,
        filtered_fasta = os.path.join(dir_reports, "mitogenome", "{sample}_consensus.fasta")
    params:
    shell("""
        set -euo pipefail

        seq_len=$(grep -v '^>' {input.fasta} | tr -d '\\n' | wc -c)
        n_count=$(grep -v '^>' {input.fasta} | tr -d '\\n' | tr 'a-z' 'A-Z' | grep -o 'N' | wc -l)

        if [ "$seq_len" -eq 0 ]; then
            echo "Empty sequence, failing QC"
            exit 1
        fi

        frac=$(awk -v n=$n_count -v l=$seq_len 'BEGIN {{print n/l}}')

        awk -v f=$frac -v m={params.max_frac} 'BEGIN {{ exit !(f <= m) }}'

        echo "PASS: N fraction = $frac"

        cp {input.fasta} {output.filtered_fasta}
        touch {output.output_fasta}
        """)