"""
Merges all the proteins commonly found in all of the samples.
So if a mitocgenome is partial and missing a gene, then that gene will not be included in the alignment and tree
"""
import os 
import glob

rule merge_proteins:
    input:
        expand(os.path.join(dir_mitos, "{sample}_mitogenome", "done.txt"), sample=sample_names)
    output:
        os.path.join(dir_out, "temp", "genes_merged.txt")
    conda:
        os.path.join(dir_env, "mafft.yaml")
    params:
        indir=os.path.join(dir_mitos),
        folder=os.path.join(dir_mitos, "mafft"),
        samles="".join(sample_names)
    shell:
        """
        mkdir -p {params.folder}
        rm -f {params.folder}/*.faa

        for gene in atp6 atp8 cob cox1 cox2 cox3 nad1 nad2 nad3 nad4 nad4l nad5 nad6
        do
            for sample in {sample_names}
            do
                f={params.indir}/${{sample}}_mitogenome")/${{sample}}_updated_result.faa_prefixed.faa.split/${{sample}}_updated_result.faa_prefixed.part_${{sample}}_${{gene}}.faa
                if [ -f "$f" ]; then
                    cat "$f" >> {params.folder}/${{gene}}.faa
                else
                    echo "File $f does not exist."
                fi
            done
        done

        touch {output}
        """

"""Aligns the merged proteins using mafft"""
rule mafft:
    input:
        expand(os.path.join(dir_out, "temp", "{sample}_merged.txt"), sample=sample_names)
    output:
        folder=os.path.join(dir_out, "temp", "aligned_done.txt")
    conda:
        os.path.join(dir_env, "mafft.yaml")
    params:
        indir=os.path.join(dir_mitos, "mafft"),
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        for f in {params.indir}/*.faa
        do
            gene=$(basename "$f" .faa)
            mafft --auto --thread {threads} "$f" > {params.indir}/${{gene}}_aligned.faa
        done
        touch {output.folder}
        """

"""Concatenates the aligned proteins into a single file for tree building"""
rule concat_alignments:
    input:
        folder=os.path.join(dir_out, "temp", "aligned_done.txt") 
    output:
        fasta=os.path.join(dir_mitos, "mafft", "concatenated_alignment.faa"),
        partitions_txt=os.path.join(dir_mitos, "mafft", "Partitions.txt"),
        partitions_nex=os.path.join(dir_mitos, "mafft", "Partitions.nex")
    params:
        script="workflow/scripts/concat_alignments.py",
        indir=os.path.join(dir_mitos, "mafft"), # folder with *_aligned.faa,
        tmp_dir=os.path.join(dir_out, "temp")
    conda:
        os.path.join(dir_env, "seqkit.yaml")  # make sure Biopython is installed
    shell:
        """
        ./{params.script} --tmp-dir {params.tmp_dir} \
            --output-fasta {output.fasta} \
            --partition-txt {output.partitions_txt} \
            --partition-nex {output.partitions_nex}
        """