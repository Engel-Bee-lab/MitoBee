"""
Rules to include host contamination statistics in final QC report
"""
rule generate_stats:
    input:
        stats=os.path.join(dir_hostcleaned,"{sample}_bamstats.txt"),
        json=os.path.join(dir_fastp,"{sample}.stats.json")
    output:
        report=os.path.join(dir_reports, "QC", "{sample}_host_qc_report.txt")
    localrule: True
    shell:
        """
        total_reads="$(
            awk '/"before_filtering"[[:space:]]*:/ {{bf=1}}
                bf && /"total_reads"[[:space:]]*:/ {{gsub(/[^0-9]/,""); print; exit}}' {input.json})"
        QC=$(cat {input.stats} | grep 'in total' | awk '{{print $1}}')
        mapped_reads=$(cat {input.stats} | grep 'mapped (' | awk '{{print $1}}' | head -n 1 )
        percent_mapped=$(cat {input.stats} | grep 'mapped (' | awk '{{print $5}}' | tr -d '()%' | head -n 1)
        
        echo "Sample: {wildcards.sample}" > {output.report}
        echo "Total Reads: $total_reads" >> {output.report}
        echo "QC Reads: $QC" >> {output.report}
        echo "Mapped Reads: $mapped_reads" >> {output.report}
        echo "Percentage of Mapped Reads: $percent_mapped%" >> {output.report}
        """

from glob import glob
sample_names = config["sample_names"]
print (sample_names)

rule final_qc_report:
    input:
        expand(os.path.join(dir_reports, "QC", "{sample}_host_qc_report.txt"), sample=sample_names)
    params:
        sample=" ".join(sample_names)
    output:
        os.path.join(dir_reports, "final_host_qc_summary.txt")
    localrule: True
    shell:
        """
        set -u  # keep undefined var checking, but don't fail on grep misses
        set -o pipefail

        echo -e "sample\ttotal_reads\tQC_reads\tmapped_reads\tpercent_mapped" > {output}

        for sample in {params.sample}; do
            report="{dir_reports}/QC/${{sample}}_host_qc_report.txt"
            total_reads=$(grep -m 1 "Total Reads:" "$report" | cut -d':' -f2 | tr -d '[:space:]' || true)
            QC_reads=$(grep -m 1 "QC Reads:" "$report" | cut -d':' -f2 | tr -d '[:space:]' || true)
            mapped=$(grep -m 1 "Mapped Reads:" "$report" | cut -d':' -f2 | tr -d '[:space:]' || true)
            percent=$(grep -m 1 "Percentage of Mapped Reads:" "$report" | cut -d':' -f2 | tr -d '%[:space:]' || true)

            # Replace empty strings with NA
            [ -n "$total_reads" ] || total_reads="NA"
            [ -n "$QC_reads" ]   || QC_reads="NA"
            [ -n "$mapped" ] || mapped="NA"
            [ -n "$percent" ] || percent="NA"

            echo -e "${{sample}}\t${{total_reads}}\t${{QC_reads}}\t${{mapped}}\t${{percent}}" >> {output}
        done
    """
