"""
Rules for host contamination removal 
"""

rule host_genes_index:
    input:
        host=config['args']['host_seq']
    output:
        index=os.path.join(dir_hostcleaned, "gene_index.1.bt2")
    params:
        host_group= os.path.join(dir_hostcleaned, "ref_genes.fasta"),
        dirs=os.path.join(dir_hostcleaned)
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
        bowtie2-build {params.host_group} {params.dirs}/gene_index
        """

rule host_gene_mapping:
    input:
        r1 = os.path.join(dir_fastp,"{sample}_R1.fastq.gz"),
        r2 = os.path.join(dir_fastp,"{sample}_R2.fastq.gz"),
        index=os.path.join(dir_hostcleaned, "gene_index")
    output:
        all_bam=os.path.join(dir_hostcleaned,"{sample}_temp.bam"),
        mapped_bam=os.path.join(dir_hostcleaned, "mitogenes", "{sample}_mapped.bam")
    params:
        bam=os.path.join(dir_hostcleaned, "gene_index")
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
        bowtie2 --very-sensitive-local -L 15 -N 1 -p {threads} \
            -x {params.bam} -1 {input.r1} -2 {input.r2} | \
            samtools sort -@ {threads} -o {output.all_bam} -

        samtools view -b -F 4 -@ {threads} -o {output.mapped_bam} {input.all_bam}
        """

from glob import glob

rule bam_mitogenes_sort:
    input:
        mapped_bam = os.path.join(dir_hostcleaned, "mitogenes", "{sample}_mapped.bam"),
    output:
        sort_bam=os.path.join(dir_hostcleaned, "mitogenes", "{sample}_mapped.sorted.bam"),
        depth=os.path.join(dir_hostcleaned, "mitogenes", "{sample}_mitogenome_snps_depth.txt")
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

rule host_mitogene_snps:
    input:
        sort_bam = os.path.join(dir_hostcleaned, "mitogenes", "{sample}_mapped.sorted.bam"),
        host= config['args']['host_seq']
    output:
        vcf = os.path.join(dir_hostcleaned, "mitogenes", "{sample}_mitogenome_snps.vcf.gz"),
        filterred_vcf = os.path.join(dir_hostcleaned, "mitogenes", "{sample}_mitogenome_snps.filtered.vcf.gz")
    params:
        prefix = os.path.join(dir_hostcleaned, "mitogenes", "{sample}_mitogene"),
        stats=os.path.join(dir_hostcleaned, "mitogenes", "{sample}_mitogene_snps.stats.txt"),
        min_depth = config['run']['min_depth'],
        qual= config['run']['qual_threshold']
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
            bcftools filter -i 'QUAL>{params.qual} && DP>{params.min_depth}' {output.vcf} -Oz -o {output.filterred_vcf}
            bcftools index {output.filterred_vcf}

            #getting the stats
            bcftools stats {output.filterred_vcf} > {params.stats}
        fi
        """

rule low_coverage_mask:
    input:
        sort_bam = os.path.join(dir_hostcleaned, "mitogenes", "{sample}_mapped.sorted.bam")
    output:
        bed = os.path.join(dir_hostcleaned, "mitogenes", "{sample}_low_cov.bed")
    conda:
        os.path.join(dir_env, "bowtie2.yaml") 
    params:
        min_depth = config['run']['min_depth'] #mask regions with less than 10x coverage
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
            | awk -v min={params.min_depth} '$3 < min {{print $1"\t"$2-1"\t"$2}}' > {output.bed}
        """

rule consensus:
    input:
        merged_vcf = os.path.join(dir_hostcleaned, "mitogenes", "{sample}_mitogenome_snps.filtered.vcf.gz"),
        host= config['args']['host_seq'],
        low_cov_bed = os.path.join(dir_hostcleaned, "mitogenes", "{sample}_low_cov.bed")
    output:
        consensus_fasta = os.path.join(dir_hostcleaned, "mitogenes", "{sample}_consensus.fasta")
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

rule qc_consensus_genes:
    input:
        fasta = os.path.join(dir_hostcleaned, "mitogenes", "{sample}_consensus.fasta")
    output:
        output_fasta = os.path.join(dir_out, "temp", "{sample}_done.txt"),
    params:
        max_frac = config['run']['max_frac'], #mask sequences with more than 5% Ns
        filtered_fasta = os.path.join(dir_reports, "mitogenes", "{sample}_consensus.fasta"),
        dirs=os.path.join(dir_reports, "mitogenes")
    shell:
        """
        set -euo pipefail
        mkdir -p {params.dirs}
        
        seq_len=$(grep -v '^>' {input.fasta} | tr -d '\\n' | wc -c)
        n_count=$(grep -v '^>' {input.fasta} | tr -d '\\n' | tr 'a-z' 'A-Z' | grep -o 'N' | wc -l)

        if [ "$seq_len" -eq 0 ]; then
            echo "Empty sequence, failing QC"
        fi

        frac=$(awk -v n=$n_count -v l=$seq_len 'BEGIN {{print n/l}}')

        if awk -v f=$frac -v m={params.max_frac} 'BEGIN {{exit !(f <= m)}}'; then
            echo "PASS: N fraction = $frac"
            cp {input.fasta} {params.filtered_fasta}
        else
            echo "FAIL: N fraction = $frac"
        fi
        touch {output.output_fasta}
        """

"""
Snakemake rules for generating mitogenome reports
"""
from glob import glob

rule mitogene_summary:
    """
    Summarize mitogenome consensus FASTA:
    filename, header, length, GC content
    """
    input:
        fasta=os.path.join(dir_hostcleaned, "mitogenes", "{sample}_consensus.fasta")
    output:
        summary=os.path.join(dir_hostcleaned, "mitogenes", "{sample}_consensus.summary.tsv")
    params:
        sample="{sample}",
        max_frac=config['run']['max_frac']
    localrule: True
    shell:
        r"""
        awk -v fname="{params.sample}_consensus.fasta" -v max_frac={params.max_frac} '
            BEGIN {{ seq="" }}
            /^>/ {{
                header = substr($0, 2)
                next
            }}
            {{
                seq = seq $0
            }}
            END {{
                len = length(seq)
                seq_upper = toupper(seq)

                gc = gsub(/[GC]/, "", seq_upper)
                n_count = gsub(/N/, "", seq_upper)
                
                gc_pct = (len > 0) ? (gc / len * 100) : 0
                n_frac = (len > 0) ? (n_count / len) : 0

                status = (n_frac <= max_frac) ? "PASS" : "FAIL"

                printf "%s\t%s\t%d\t%.2f\t%d\t%.4f\t%s\n", \
                    fname, header, len, gc_pct, n_count, n_frac, status
            }}
        ' {input.fasta} > {output.summary}
        """

rule mitogene_reports_aggregate:
    input:
        summaries = expand(os.path.join(dir_hostcleaned, "mitogenes", "{sample}_consensus.summary.tsv"), sample=samples)
    output:
        aggregate = os.path.join(dir_reports, "mitogene_consensus_summary.tsv")
    localrule: True
    shell:
        """
        echo -e "filename\theader\tlength\tGC_content\tN_count\tN_fraction\tQC_status" > {output.aggregate}
        cat {input.summaries} >> {output.aggregate}
        """