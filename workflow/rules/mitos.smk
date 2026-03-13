rule install_database:
    output:
        os.path.join(dir_out, "database", "mitos_db" "mitos_downlaoded.txt")
    params:
        url ="https://zenodo.org/api/records/4284483/files-archive",
        out = os.path.join(dir_out, "database", "mitos_db.zip"),
        decom=os.path.join(dir_out, "database", "mitos_db")
    shell:
        """
        set -euo pipefail
        wget -O {params.out} {params.url}

        unzip {params.out}
        rm -rf {params.decom}
        mkdir -p {params.decom}
        for f in refseq*; do
            tar -xvjf "$f"
        done
        mv refseq* {params.decom}/.
        touch {output}
        """

rule run_mitos:
    input:
        fasta= lambda wc: os.path.join(input_dir, f"{wc.sample}.{extn}"),
        db=os.path.join(dir_out, "database", "mitos_db" "mitos_downlaoded.txt")
    output:
        os.path.join(dir_mitos, "{sample}_mitogenome", "{sample}_result.faa")
    params:
        outdir=os.path.join(dir_mitos, "{sample}_mitogenome"),
        genetic_code=2,
        database=os.path.join(dir_out, "database", "mitos_db"),
        specific="refseq89m",
        sample="{sample}"
    conda:
        os.path.join(dir_env, "mitos.yaml")
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        set -euo pipefail
        rm -rf {params.outdir}
        mkdir {params.outdir}
        runmitos -i {input.fasta} -o {params.outdir} -c {params.genetic_code} -r {params.database}/{params.specific}
        for file in {params.outdir}/*; do
            base=$(basename "$file")
            mv "$file" "{params.outdir}/{params.sample}_${{base}}"
        done
        """

rule separate_faa:
    input:
        faa=os.path.join(dir_mitos, "{sample}_mitogenome", "{sample}_result.faa")
    output:
        os.path.join(dir_mitos, "{sample}_mitogenome", "done.txt")
    params:
        outdir=os.path.join(dir_mitos, "{sample}_mitogenome"),
        updatedfaa=os.path.join(dir_mitos, "{sample}_mitogenome", "{sample}_updated_result.faa"),
    conda:
        os.path.join(dir_env, "seqkit.yaml")
    shell:
        """
        set -euo pipefail
        seqkit replace -p '.*; ' -r '' {input.faa} -o {params.updatedfaa}
        seqkit split -i {params.updatedfaa}
        touch {output}
        """
