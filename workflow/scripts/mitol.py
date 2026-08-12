#!/usr/bin/env python3
"""
Script shared from Dr. Sunagawa, on assembling mitochondrial genomes from metagenomic data. 
Being adopted to this workflow, to use the assembly function while retaining the current workflows quality control and mapping options.
"""
import argparse
import concurrent.futures
import logging
import re
import shlex
import subprocess
import sys
from pathlib import Path
import shutil

VERSION = "5.3"
DRY_RUN = False  # set in main() before any mode function runs


# ---------------- CORE HELPERS ---------------- #

def run(cmd, **kwargs):
    cmd = [str(x) for x in cmd]
    logging.debug("CMD: %s", " ".join(cmd))
    if not DRY_RUN:
        subprocess.run(cmd, check=True, **kwargs)


def shell_run(cmd):
    logging.debug("CMD: %s", cmd)
    if not DRY_RUN:
        subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True)


def require(*binaries):
    for b in binaries:
        if shutil.which(b) is None:
            logging.error("Required executable '%s' not found in PATH.", b)
            sys.exit(1)


def check_file(path, label=""):
    if not Path(path).is_file():
        logging.error("Missing %s: %s", label or "file", path)
        sys.exit(1)


def validate_sample_name(sample):
    if not re.fullmatch(r"[\w.\-]+", sample):
        logging.error(
            "Invalid sample name '%s'. Only alphanumerics, dots, underscores, and hyphens are allowed.",
            sample,
        )
        sys.exit(1)


def nonempty(path):
    p = Path(path)
    return p.exists() and p.stat().st_size > 0


def safe_copy(src, dst):
    logging.debug("COPY: %s -> %s", src, dst)
    if not DRY_RUN:
        shutil.copy(src, dst)


def shell_run_log(cmd, log_file):
    """
    Run a shell pipeline with pipefail and redirect ALL stderr (from every
    command in the pipeline) to *log_file* in append mode.
    Wrapping in { ...; } ensures the redirection applies to the whole group.
    """
    full = f"{{ {cmd}; }} 2>> {shlex.quote(str(log_file))}"
    logging.debug("CMD: %s", full)
    if not DRY_RUN:
        subprocess.run(["bash", "-o", "pipefail", "-c", full], check=True)


def has_reads_gz(path):
    """
    Return True if a gzipped FASTQ file contains at least one read.
    gzip -c of an empty stream produces a ~20-byte envelope that passes a
    simple file-size check; this function reads the first byte through gzip
    to distinguish a real file from an empty envelope.
    """
    import gzip
    if DRY_RUN:
        return True
    try:
        with gzip.open(str(path), "rb") as fh:
            return bool(fh.read(1))
    except Exception:
        return False


def pident_awk(cutoff):
    """
    Return an awk command that filters SAM alignment lines (no header) by
    pairwise percent identity:

        identity = (alignment_length - NM) / alignment_length * 100

    alignment_length is derived from the CIGAR string (sum of M, I, D
    operations), which excludes soft/hard-clipped bases from the denominator.
    Using raw query length (length($10)) instead would include clipped bases
    and artificially lower the identity, over-filtering real mito reads.

    Pair this with `samtools view -h -F 2308` so that (a) the SAM header is
    preserved for downstream `samtools fastq`, and (b) secondary alignments
    (which carry SEQ=* and would corrupt the calculation) are already gone.
    """
    return (
        f"awk '{{"
        f"if($1~/^@/){{print;next}};"  # pass SAM header lines through unchanged
        f"alen=0; c=$6; n=\"\";"
        f" for(k=1;k<=length(c);k++){{"
        f"ch=substr(c,k,1);"
        f" if(ch~/[0-9]/) n=n ch;"
        f" else{{if(ch~/[MID]/) alen+=n+0; n=\"\"}}"
        f"}};"
        f" nm=0;"
        f" for(i=12;i<=NF;i++) if($i~/^NM:i:/) {{nm=substr($i,6)+0; break}};"
        f" if(alen>0 && (alen-nm)/alen*100>={cutoff}) print"
        f"}}'"
    )


def find_fastq(directory, sample, read_num):
    """
    Locate a FASTQ file for *sample* and *read_num* (1 or 2) under *directory*.
    Tries dot and underscore separators, and both .fq / .fastq with or without .gz.
    Returns a resolved Path, or None if nothing matches.
    """
    for sep in (".", "_"):
        for ext in ("fq.gz", "fastq.gz", "fq", "fastq"):
            p = directory / f"{sample}{sep}{read_num}.{ext}"
            if p.exists():
                return p.resolve()
    return None


def find_bbtools_resource(filename):
    """
    Locate a file in the BBTools bundled resources directory.
    BBTools ships adapters.fa and phix174_ill.ref.fa.gz alongside bbduk.sh,
    so as long as bbduk.sh is in PATH the files are always present — no
    separate resources/ directory needed.
    Returns a resolved Path, or None if not found.
    """
    bbduk = shutil.which("bbduk.sh")
    if bbduk is None:
        return None
    bbduk_dir = Path(bbduk).resolve().parent
    for candidate in [
        bbduk_dir / "resources" / filename,   # standard BBTools layout
        bbduk_dir / filename,                  # flat layout (some conda installs)
        bbduk_dir.parent / "resources" / filename,  # wrapper-script layout
    ]:
        if candidate.is_file():
            return candidate
    return None


# ---------------- QC ---------------- #

def qc_mode(args):
    require("bbduk.sh")
    validate_sample_name(args.sample)

    sample    = args.sample
    reads_dir = Path(args.reads_dir).resolve()
    out_dir   = Path("qc.out") / sample
    out_dir.mkdir(parents=True, exist_ok=True)

    log_file = out_dir / args.log
    # Truncate log on a fresh run so re-runs don't silently accumulate output
    if not args.append_log and log_file.exists() and not DRY_RUN:
        log_file.unlink()

    fq1 = find_fastq(reads_dir, sample, 1)
    fq2 = find_fastq(reads_dir, sample, 2)
    for fq, label in [(fq1, "R1"), (fq2, "R2")]:
        if fq is None:
            logging.error(
                "Could not find %s FASTQ for sample '%s' in %s", label, sample, reads_dir
            )
            sys.exit(1)

    adapters_ref = Path(args.adapters).resolve() if args.adapters else find_bbtools_resource("adapters.fa")
    phix_ref     = Path(args.phix).resolve()     if args.phix     else find_bbtools_resource("phix174_ill.ref.fa.gz")

    if adapters_ref is None:
        logging.error(
            "Could not locate adapters.fa in the BBTools installation. "
            "Pass --adapters to specify the path explicitly."
        )
        sys.exit(1)
    if phix_ref is None:
        logging.error(
            "Could not locate phix174_ill.ref.fa.gz in the BBTools installation. "
            "Pass --phix to specify the path explicitly."
        )
        sys.exit(1)

    check_file(adapters_ref, "adapters reference")
    check_file(phix_ref,     "PhiX reference")
    logging.debug("Adapters ref: %s", adapters_ref)
    logging.debug("PhiX ref:     %s", phix_ref)

    clean_r1 = out_dir / f"{sample}.clean.1.fq.gz"
    clean_r2 = out_dir / f"{sample}.clean.2.fq.gz"
    clean_s  = out_dir / f"{sample}.clean.s.fq.gz"
    clean_m  = out_dir / f"{sample}.clean.m.fq.gz"

    mem = args.mem
    q   = shlex.quote
    cmd = (
        f"bbduk.sh -Xmx{mem} usejni=t in={q(str(fq1))} in2={q(str(fq2))} "
        f"out=stdout.fq "
        f"outm={q(str(out_dir / f'{sample}_adapter_matched.fq.gz'))} "
        f"outs={q(str(out_dir / f'{sample}_adapter_s.fq.gz'))} "
        f"refstats={q(str(out_dir / f'{sample}.adapter_trim.stats'))} "
        f"statscolumns=5 overwrite=t ref={q(str(adapters_ref))} 2>> {q(str(log_file))} | "
        f"bbduk.sh -Xmx{mem} usejni=t interleaved=true overwrite=t "
        f"in=stdin.fq out=stdout.fq "
        f"outm={q(str(out_dir / f'{sample}_phix_matched.fq.gz'))} "
        f"outs={q(str(out_dir / f'{sample}_phix_s.fq.gz'))} "
        f"ref={q(str(phix_ref))} k=31 hdist=1 "
        f"refstats={q(str(out_dir / f'{sample}_phix.stats'))} "
        f"statscolumns=5 2>> {q(str(log_file))} | "
        f"bbduk.sh -Xmx{mem} usejni=t overwrite=t interleaved=true "
        f"in=stdin.fq out1={q(str(clean_r1))} out2={q(str(clean_r2))} "
        f"outm={q(str(clean_m))} outs={q(str(clean_s))} minlength=45 "
        f"qtrim=rl maq=20 maxns=1 "
        f"stats={q(str(out_dir / f'{sample}_qc.stats'))} "
        f"statscolumns=5 trimq=14 2>> {q(str(log_file))}"
    )
    shell_run(cmd)


# ---------------- MITO ASSEMBLY ---------------- #

def assemble_mt_mode(args):
    # repair.sh and sam2fastq.pl are no longer required; samtools and seqkit replace them
    require("bwa", "samtools", "seqkit", "spades.py")
    validate_sample_name(args.sample_id)

    q         = shlex.quote
    sample_id = args.sample_id
    mito_db   = Path(args.mitodb).resolve()

    base_output_dir = (
        Path(args.output_dir).resolve()
        if args.output_dir
        else Path("assembly.out").resolve() / sample_id
    )
    base_output_dir.mkdir(parents=True, exist_ok=True)

    read1  = Path(args.read1).resolve()
    read2  = Path(args.read2).resolve()
    single = Path(args.single).resolve() if args.single else None
    merged = Path(args.merged).resolve() if args.merged else None

    check_file(mito_db, "mitochondrial reference")
    check_file(read1, "forward reads")
    check_file(read2, "reverse reads")
    if single:
        check_file(single, "single reads")
    if merged:
        check_file(merged, "merged reads")

    tmp_dir  = base_output_dir / f"tmp_{sample_id}"
    run_log  = base_output_dir / f"{sample_id}.assemble_mt.log"

    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    if not DRY_RUN:
        run_log.unlink(missing_ok=True)   # fresh log each run

    logging.info("Run log: %s", run_log)

    # Pre-declare so they're always in scope for the summary block after try/finally
    synced_r1    = base_output_dir / f"{sample_id}.synced.1.fq.gz"
    synced_r2    = base_output_dir / f"{sample_id}.synced.2.fq.gz"
    final_single = None
    final_merged = None
    final_fasta  = base_output_dir / f"{sample_id}.fasta"
    final_mode   = None
    spades_out   = None

    try:
        def slog(cmd):
            """shell_run_log bound to this run's log file."""
            shell_run_log(cmd, run_log)

        def map_to_names(reads, out_names, n_threads):
            """
            Map reads single-end, collect QNAMEs of alignments passing the
            percent-identity cutoff.  Strips legacy /1 /2 suffixes so names
            are compatible with both R1 and R2 originals.
            Assumes modern Illumina format (R1 and R2 share identical names).
            """
            slog(
                f"bwa mem -t {n_threads} {q(str(mito_db))} {q(str(reads))} "
                f"| samtools view -h -F 2308 "
                f"| {pident_awk(args.cutoff)} "
                f"| cut -f1 | sed 's|/[12]$||' | sort -u > {q(str(out_names))}"
            )

        def count_lines(path):
            """Return the number of lines in a text file (= mapped read count)."""
            if DRY_RUN or not Path(path).exists():
                return 0
            return sum(1 for _ in open(path))

        def run_spades(mode, s_r1, s_r2, fs=None, fm=None):
            out           = base_output_dir / f"spades_{mode}_out"
            spades_log    = base_output_dir / f"spades_{mode}.log"

            if out.exists():
                logging.warning("Removing existing SPAdes output directory: %s", out)
                shutil.rmtree(out)

            cmd = ["spades.py", f"--{mode}", "-1", s_r1, "-2", s_r2, "-t", args.threads]
            if fs and has_reads_gz(fs):
                cmd.extend(["-s", fs])
            if fm and has_reads_gz(fm):
                cmd.extend(["--merged", fm])
            cmd.extend(["-o", out])

            logging.info("Running SPAdes --%s  (log: %s)", mode, spades_log)
            if not DRY_RUN:
                with open(spades_log, "w") as log:
                    subprocess.run(
                        [str(x) for x in cmd],
                        stdout=log, stderr=subprocess.STDOUT, check=True,
                    )
            return out

        bwa_index_files = [f"{mito_db}.{ext}" for ext in ["amb", "ann", "bwt", "pac", "sa"]]
        if not all(Path(f).is_file() for f in bwa_index_files):
            logging.info("Creating BWA index...")
            slog(f"bwa index {q(str(mito_db))}")
        else:
            logging.info("BWA index already exists. Skipping.")

        # Map R1 and R2 in parallel using half the threads each to avoid
        # over-subscription on the calling node.
        par_threads = max(1, args.threads // 2)
        names_r1  = tmp_dir / "names_r1.txt"
        names_r2  = tmp_dir / "names_r2.txt"
        all_names = tmp_dir / "all_names.txt"

        logging.info("Mapping R1 and R2 against mito reference (parallel, %d threads each)...", par_threads)
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            f1 = pool.submit(map_to_names, read1, names_r1, par_threads)
            f2 = pool.submit(map_to_names, read2, names_r2, par_threads)
            f1.result()
            f2.result()

        n_r1, n_r2 = count_lines(names_r1), count_lines(names_r2)
        logging.info("Mapped read names: R1=%d  R2=%d", n_r1, n_r2)

        # Paired rescue: union of R1 and R2 mapped names only.
        # Merged reads are NOT included here — merged reads and their source
        # R1/R2 represent the same fragment; cross-rescuing would duplicate
        # those fragments for the assembler.
        shell_run(
            f"cat {q(str(names_r1))} {q(str(names_r2))} "
            f"| sort -u > {q(str(all_names))}"
        )
        logging.info("Union of mapped names: %d", count_lines(all_names))

        logging.info("Extracting read pairs (including rescued mates) from original FASTQs...")
        run(["seqkit", "grep", "-f", str(all_names), str(read1), "-o", str(synced_r1)])
        run(["seqkit", "grep", "-f", str(all_names), str(read2), "-o", str(synced_r2)])

        if merged:
            # Map merged reads independently and keep only those that align;
            # no cross-rescue with the paired reads.
            final_merged = base_output_dir / f"{sample_id}.mrg.fq.gz"
            slog(
                f"bwa mem -t {args.threads} {q(str(mito_db))} {q(str(merged))} "
                f"| samtools view -h -F 2308 "
                f"| {pident_awk(args.cutoff)} "
                f"| samtools fastq -0 /dev/stdout - "
                f"| gzip -c > {q(str(final_merged))}"
            )
            if not has_reads_gz(final_merged):
                logging.warning(
                    "No merged reads passed the %.0f%% identity filter — "
                    "%s will not be passed to SPAdes.", args.cutoff, final_merged
                )
                final_merged = None

        if single:
            final_single = base_output_dir / f"{sample_id}.single.fq.gz"
            slog(
                f"bwa mem -t {args.threads} {q(str(mito_db))} {q(str(single))} "
                f"| samtools view -h -F 2308 "
                f"| {pident_awk(args.cutoff)} "
                f"| samtools fastq -0 /dev/stdout - "
                f"| gzip -c > {q(str(final_single))}"
            )
            if not has_reads_gz(final_single):
                logging.warning(
                    "No single reads passed the %.0f%% identity filter — "
                    "%s will not be passed to SPAdes.", args.cutoff, final_single
                )
                final_single = None

        for mode in ["plasmid", "metaplasmid", "meta"]:
            spades_out   = run_spades(mode, synced_r1, synced_r2, final_single, final_merged)
            contigs_file = spades_out / "contigs.fasta"

            if nonempty(contigs_file):
                safe_copy(contigs_file, final_fasta)
                final_mode = mode
                logging.info("Wrote candidate assembly: %s", final_fasta)
                break

            logging.info("No non-empty contigs.fasta from --%s", mode)

    except subprocess.CalledProcessError as e:
        logging.error("Command failed: %s", e)
        sys.exit(1)

    finally:
        if tmp_dir.exists():
            shutil.rmtree(tmp_dir)

    logging.info("Done.")
    logging.info("  Synced reads:    %s, %s", synced_r1, synced_r2)
    if final_single and final_single.exists():
        logging.info("  Single reads:    %s", final_single)
    if final_merged and final_merged.exists():
        logging.info("  Merged reads:    %s", final_merged)
    if final_mode:
        logging.info("  Final mode:      --%s", final_mode)
    if spades_out:
        logging.info("  Assembly dir:    %s", spades_out)
    if final_fasta.exists():
        logging.info("  Candidate FASTA: %s", final_fasta)
    else:
        logging.warning("No candidate FASTA was produced.")


# ---------------- EXTRACT BEST MITO CONTIG ---------------- #

def extract_mode(args):
    require("minimap2", "seqkit")

    sample = args.sample
    ref    = Path(args.ref).resolve()
    check_file(ref, "mitochondrial reference")

    if args.input:
        input_fasta = Path(args.input).resolve()
    else:
        input_dir   = Path(args.input_dir).resolve() if args.input_dir else Path(f"{sample}.mt").resolve()
        input_fasta = input_dir / f"{sample}.fasta"

    check_file(input_fasta, "input candidate FASTA")

    output_fasta = (
        Path(args.output).resolve()
        if args.output
        else input_fasta.parent / f"{sample}.mt.fasta"
    )

    if output_fasta.exists() and not args.force:
        logging.error("Output already exists: %s  (use --force to overwrite)", output_fasta)
        sys.exit(1)

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    mm2_log = output_fasta.parent / f"{sample}.minimap2.log"

    logging.info("Sample:        %s", sample)
    logging.info("Input FASTA:   %s", input_fasta)
    logging.info("Reference:     %s", ref)
    logging.info("Output FASTA:  %s", output_fasta)
    logging.info("minimap2 log:  %s", mm2_log)

    result_stdout = ""
    if not DRY_RUN:
        try:
            with open(mm2_log, "w") as mlog:
                result = subprocess.run(
                    ["minimap2", "-t", str(args.threads), "-x", args.preset,
                     str(input_fasta), str(ref)],
                    check=True, text=True,
                    stdout=subprocess.PIPE, stderr=mlog,
                )
            result_stdout = result.stdout
        except subprocess.CalledProcessError as e:
            logging.error("minimap2 failed: %s", e)
            sys.exit(1)

    best = None
    for line in result_stdout.splitlines():
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 12:
            continue
        tags = fields[12:]
        if any(tag == "tp:A:S" for tag in tags):
            continue

        s1 = None
        for tag in tags:
            if tag.startswith("s1:i:"):
                try:
                    s1 = int(tag.split(":")[-1])
                except ValueError:
                    pass

        if s1 is None:
            continue

        record = {"s1": s1, "query": fields[0], "strand": fields[4], "target": fields[5]}
        if best is None or record["s1"] > best["s1"]:
            best = record

    if best is None:
        logging.warning("No primary minimap2 hit found for %s", sample)
        sys.exit(0)

    contig = best["target"]
    strand = best["strand"]

    logging.info("Best contig:   %s", contig)
    logging.info("Strand:        %s", strand)
    logging.info("s1 score:      %d", best["s1"])

    if args.min_score is not None and best["s1"] < args.min_score:
        logging.warning(
            "Best s1 score %d is below --min-score %d. Skipping.",
            best["s1"], args.min_score,
        )
        sys.exit(0)

    # seqkit uses Go/RE2 (no POSIX classes); match ID field only so $ suffices
    pattern    = rf"^{re.escape(contig)}$"
    seqkit_cmd = ["seqkit", "grep", "-r", "-p", pattern, str(input_fasta)]
    seq_cmd    = (
        ["seqkit", "seq", "-r", "-p", "-w", "0"]
        if args.orient_to_ref and strand == "-"
        else ["seqkit", "seq", "-w", "0"]
    )

    try:
        p1 = subprocess.Popen(seqkit_cmd, stdout=subprocess.PIPE, text=True)
        p2 = subprocess.Popen(seq_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, text=True)

        if p1.stdout is None:
            logging.error("Failed to open pipe to seqkit grep")
            sys.exit(1)

        p1.stdout.close()
        seq_output, _ = p2.communicate()
        p1_ret = p1.wait()

        if p1_ret != 0 or p2.returncode != 0:
            logging.error("seqkit extraction failed")
            sys.exit(1)
    except Exception as e:
        logging.error("Failed during seqkit extraction: %s", e)
        sys.exit(1)

    if not seq_output.strip():
        logging.error("seqkit extracted no sequence for contig: %s", contig)
        sys.exit(1)

    lines = seq_output.splitlines()
    if not lines or not lines[0].startswith(">"):
        logging.error("Extracted sequence does not look like FASTA")
        sys.exit(1)

    lines[0] = f">{sample} {lines[0][1:]}"

    if not DRY_RUN:
        with open(output_fasta, "w") as out:
            out.write("\n".join(lines) + "\n")

    logging.info("Wrote: %s", output_fasta)


# ---------------- DENOVO ---------------- #

def assemble_denovo_mode(args):
    require("spades.py")
    validate_sample_name(args.sample)

    sample  = args.sample
    qc_base = Path(args.qc_dir).resolve() / sample

    read1  = qc_base / f"{sample}.clean.1.fq.gz"
    read2  = qc_base / f"{sample}.clean.2.fq.gz"
    single = qc_base / f"{sample}.clean.s.fq.gz"
    merged = qc_base / f"{sample}.clean.m.fq.gz"

    check_file(read1, "QC forward reads")
    check_file(read2, "QC reverse reads")

    base_output_dir = (
        Path(args.output_dir).resolve()
        if args.output_dir
        else Path("metaspades_denovo_out").resolve() / sample
    )
    base_output_dir.mkdir(parents=True, exist_ok=True)

    spades_out = base_output_dir / "spades_meta_out"
    log_file   = base_output_dir / "spades_meta.log"

    cmd = [
        "spades.py", "--meta",
        "-1", str(read1), "-2", str(read2),
        "-t", str(args.threads),
        "-o", str(spades_out),
    ]
    if nonempty(merged):
        cmd.extend(["--merged", str(merged)])
    if nonempty(single):
        cmd.extend(["-s", str(single)])

    if spades_out.exists():
        logging.warning("Removing existing SPAdes output directory: %s", spades_out)
        shutil.rmtree(spades_out)

    logging.info("Running SPAdes --meta  (log: %s)", log_file)
    try:
        if not DRY_RUN:
            with open(log_file, "w") as log:
                subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
    except subprocess.CalledProcessError as e:
        logging.error("SPAdes failed: %s", e)
        sys.exit(1)

    logging.info("Done.  Assembly dir: %s", spades_out)


# ---------------- TRIM ---------------- #

def _trim_redundant_overlap(seq, min_overlap, max_overlap):
    max_possible = min(max_overlap, len(seq) // 2)
    for i in range(max_possible, min_overlap - 1, -1):
        if seq[:i] == seq[-i:]:
            return seq[:-i], i
    return seq, 0


def trim_mode(args):
    try:
        from Bio import SeqIO
        from Bio.Seq import Seq
    except ImportError:
        logging.error("Biopython is required for 'trim'.  Install with: pip install biopython")
        sys.exit(1)

    sample   = args.sample
    mt_dir   = Path(f"{sample}.mt")
    in_path  = Path(args.input)  if args.input  else mt_dir / f"{sample}.mt.fasta"
    out_path = Path(args.output) if args.output else mt_dir / f"{sample}.mt.trimmed.fasta"

    logging.info("Input : %s", in_path)
    logging.info("Output: %s", out_path)

    total   = 0
    trimmed = 0
    records = []

    with open(in_path) as fh:
        for record in SeqIO.parse(fh, "fasta"):
            total += 1
            original_len = len(record.seq)
            seq_str, cut_len = _trim_redundant_overlap(
                str(record.seq), args.min_overlap, args.max_overlap
            )
            new_len = len(seq_str)

            if cut_len > 0:
                logging.info("%s: trimmed %d bp (%d -> %d)", record.id, cut_len, original_len, new_len)
                trimmed += 1

            record.seq = Seq(seq_str)
            # Re-derive id and description independently to avoid BioPython
            # writing a duplicate id prefix in the FASTA header.
            updated_full = re.sub(r"length_\d+", f"length_{new_len}", record.description)
            record.id          = updated_full.split()[0]
            record.description = " ".join(updated_full.split()[1:])
            records.append(record)

    if not DRY_RUN:
        with open(out_path, "w") as out:
            SeqIO.write(records, out, "fasta")

    logging.info("Processed %d sequences; trimmed %d overhangs.", total, trimmed)


# ---------------- ROTATE ---------------- #

def rotate_mode(args):
    try:
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
    except ImportError:
        logging.error("Biopython is required for 'rotate'.  Install with: pip install biopython")
        sys.exit(1)

    def _rotate(seq_str, start):
        return seq_str[start:] + seq_str[:start]

    def _find_and_orient(seq_str, regex):
        m = regex.search(seq_str)
        if m:
            return Seq(_rotate(seq_str, m.start())), "forward"
        rc = str(Seq(seq_str).reverse_complement())
        m  = regex.search(rc)
        if m:
            return Seq(_rotate(rc, m.start())), "reverse"
        return Seq(seq_str), "unmatched"

    sample   = args.sample
    mt_dir   = Path(f"{sample}.mt")
    in_path  = Path(args.input)  if args.input  else mt_dir / f"{sample}.mt.trimmed.fasta"
    out_path = Path(args.output) if args.output else mt_dir / f"{sample}.mt.trimmed.rotated.fasta"
    log_path = Path(args.log)    if args.log    else mt_dir / f"{sample}.rotation_log.tsv"

    logging.info("Input : %s", in_path)
    logging.info("Output: %s", out_path)
    logging.info("Log   : %s", log_path)

    motif_regex = re.compile(args.motif, re.IGNORECASE)
    records     = []
    log_entries = []

    for record in SeqIO.parse(in_path, "fasta"):
        logging.info("Processing %s", record.id)
        new_seq, orientation = _find_and_orient(str(record.seq), motif_regex)
        # Strip the leading id from description to avoid BioPython writing a
        # duplicate id prefix in the output FASTA header.
        desc = record.description[len(record.id):].lstrip()
        records.append(SeqRecord(new_seq, id=record.id, description=desc))
        log_entries.append((record.id, orientation))

    unmatched_ids = [rid for rid, o in log_entries if o == "unmatched"]
    if unmatched_ids:
        logging.warning(
            "%d sequence(s) had no motif match and were left unchanged: %s",
            len(unmatched_ids),
            ", ".join(unmatched_ids),
        )
        if args.strict:
            logging.error("Aborting due to unmatched sequences (--strict).")
            sys.exit(1)

    if not DRY_RUN:
        if records:
            SeqIO.write(records, out_path, "fasta")
            logging.info("Wrote %d sequences to %s", len(records), out_path)

        with open(log_path, "w") as logf:
            logf.write("record_id\torientation\n")
            for rec_id, orient in log_entries:
                logf.write(f"{rec_id}\t{orient}\n")
        logging.info("Orientation log saved to %s", log_path)


# ---------------- MAIN ---------------- #

def main():
    global DRY_RUN

    # Shared parent parser for subcommands that use --threads
    _tp = argparse.ArgumentParser(add_help=False)
    _tp.add_argument("--threads", type=int, default=32, help="CPU threads [default: 32]")

    parser = argparse.ArgumentParser(description="miTOL toolkit")
    parser.add_argument("--version",  action="version", version=f"miTOL {VERSION}")
    parser.add_argument("--dry-run",  action="store_true",
                        help="Print commands without executing them")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Show commands as they run (DEBUG output)")

    sub = parser.add_subparsers(dest="command")

    # qc
    p = sub.add_parser("qc", help="Adapter/PhiX trimming and quality filtering")
    p.add_argument("sample")
    p.add_argument("--reads-dir",  default=".",      help="Directory containing raw FASTQ files [default: .]")
    p.add_argument("--log",        default="qc.log", help="Log filename inside output directory [default: qc.log]")
    p.add_argument("--mem",        default="1G",     help="Memory per BBDuk step, e.g. 4G [default: 1G]")
    p.add_argument("--adapters",   default=None,
                   help="Adapters FASTA [default: auto-detected from BBTools installation]")
    p.add_argument("--phix",       default=None,
                   help="PhiX reference FASTA [default: auto-detected from BBTools installation]")
    p.add_argument("--append-log", action="store_true",
                   help="Append to existing log instead of truncating")
    p.set_defaults(func=qc_mode)

    # assemble_mt
    p = sub.add_parser("assemble_mt", parents=[_tp],
                       help="Map reads to a mito reference and assemble with SPAdes")
    p.add_argument("sample_id")
    p.add_argument("mitodb", help="Mitochondrial reference FASTA (BWA-indexed if needed)")
    p.add_argument("-1", dest="read1", required=True, help="Forward reads (R1)")
    p.add_argument("-2", dest="read2", required=True, help="Reverse reads (R2)")
    p.add_argument("-m", dest="merged", help="Merged reads")
    p.add_argument("-s", dest="single", help="Single-end reads")
    p.add_argument("--output-dir", help="Output directory [default: assembly.out/{sample_id}]")
    p.add_argument("--cutoff", type=int, default=90,
                   help="Minimum percent identity to keep a mapped read [default: 90]")
    p.set_defaults(func=assemble_mt_mode)

    # extract
    p = sub.add_parser("extract", parents=[_tp],
                       help="Extract best mito contig from a candidate assembly")
    p.add_argument("sample", help="Sample ID")
    p.add_argument("ref",    help="Reference mitochondrial FASTA")
    p.add_argument("-i", "--input",     help="Input candidate FASTA [default: ./SAMPLE.mt/SAMPLE.fasta]")
    p.add_argument("--input-dir",       help="Input directory [default: ./SAMPLE.mt]")
    p.add_argument("-o", "--output",    help="Output FASTA [default: input dir/SAMPLE.mt.fasta]")
    p.add_argument("--preset",          default="asm20", help="minimap2 preset [default: asm20]")
    p.add_argument("--min-score",       type=int, default=None,
                   help="Minimum s1:i score required to write output")
    p.add_argument("--orient-to-ref",   action="store_true", default=True,
                   help="Reverse-complement if best hit is on minus strand [default: on]")
    p.add_argument("--no-orient-to-ref", dest="orient_to_ref", action="store_false",
                   help="Disable reverse-complement of minus-strand hits")
    p.add_argument("--force",           action="store_true",
                   help="Overwrite existing output FASTA")
    p.set_defaults(func=extract_mode)

    # assemble_denovo
    p = sub.add_parser("assemble_denovo", parents=[_tp],
                       help="Full metaSPAdes de-novo assembly from QC reads")
    p.add_argument("sample")
    p.add_argument("--qc-dir",     default="qc.out",
                   help="Root QC output directory [default: qc.out]")
    p.add_argument("--output-dir", help="Output directory [default: metaspades_denovo_out/{sample}]")
    p.set_defaults(func=assemble_denovo_mode)

    # trim
    p = sub.add_parser("trim", help="Trim repeated prefix/suffix overhangs from circular sequences")
    p.add_argument("sample", help="Sample name")
    p.add_argument("-i", "--input",
                   help="Input FASTA [default: SAMPLE.mt/SAMPLE.mt.fasta]")
    p.add_argument("-o", "--output",
                   help="Output FASTA [default: SAMPLE.mt/SAMPLE.mt.trimmed.fasta]")
    p.add_argument("--min-overlap", type=int, default=10,
                   help="Minimum overlap length [default: 10]")
    p.add_argument("--max-overlap", type=int, default=1000,
                   help="Maximum overlap length [default: 1000]")
    p.set_defaults(func=trim_mode)

    # rotate
    p = sub.add_parser("rotate", help="Reorient circular sequences to start at a regex motif")
    p.add_argument("sample", help="Sample name")
    p.add_argument("-i", "--input",
                   help="Input FASTA [default: SAMPLE.mt/SAMPLE.mt.trimmed.fasta]")
    p.add_argument("-o", "--output",
                   help="Output FASTA [default: SAMPLE.mt/SAMPLE.mt.trimmed.rotated.fasta]")
    p.add_argument("--log",
                   help="TSV orientation log [default: SAMPLE.mt/SAMPLE.rotation_log.tsv]")
    p.add_argument("--motif",  default=r"CCAGC[AC]GACGCGGT[AG]A[ACGT]ACTTA",
                   help="Regex motif to orient to [default: built-in mitochondrial motif]")
    p.add_argument("--strict", action="store_true",
                   help="Exit non-zero if any sequence has no motif match")
    p.set_defaults(func=rotate_mode)

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    logging.basicConfig(
        format="[%(levelname)s] %(message)s",
        level=logging.DEBUG if args.verbose else logging.INFO,
    )

    DRY_RUN = args.dry_run
    if DRY_RUN:
        logging.info("Dry-run mode: commands will be printed but not executed.")

    args.func(args)


if __name__ == "__main__":
    main()
