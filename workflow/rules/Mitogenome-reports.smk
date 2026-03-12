"""
Snakemake rules for generating mitogenome reports
"""
from glob import glob

rule mitogenome_summary:
    """
    Summarize mitogenome consensus FASTA:
    filename, header, length, GC content
    """
    input:
        fasta=os.path.join(dir_hostcleaned, "mitogenome", "{sample}_consensus.fasta")
    output:
        summary=os.path.join(dir_hostcleaned, "mitogenome", "{sample}_consensus.summary.tsv")
    params:
        sample="{sample}",
        max_frac=0.333
    localrule: True
    shell:
        r"""
        awk -v fname="{params.sample}_consensus.fasta" -v max_frac={params.max_frac} '
            BEGIN { seq="" }
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

rule mitogenome_reports_aggregate:
    input:
        summaries = os.path.join(dir_hostcleaned, "mitogenome", "*_consensus.summary.tsv")
    output:
        aggregate = os.path.join(dir_reports, "mitogenome_consensus_summary.tsv")
    localrule: True
    shell:
        """
        echo -e "filename\theader\tlength\tGC_content\tN_count\tN_fraction\tQC_status" > {output.aggregate}
        cat {input.summaries} >> {output.aggregate}
        """