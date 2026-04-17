import pandas as pd

# Snakemake inputs
primary_idx = snakemake.input.primary_idx
primary_cov = snakemake.input.primary_cov
strict_idx = snakemake.input.strict_idx
output_file = snakemake.output.summary

# -----------------------------
# Load idxstats (primary)
# -----------------------------
primary = pd.read_csv(
    primary_idx,
    sep="\t",
    header=None,
    names=["reference", "length", "primary_reads", "unmapped"]
)

primary = primary[primary["reference"] != "*"]

# -----------------------------
# Load strict (MAPQ ≥ 30)
# -----------------------------
strict = pd.read_csv(
    strict_idx,
    sep="\t",
    header=None,
    names=["reference", "length2", "strict", "unmapped2"]
)

strict = strict[["reference", "strict"]]

# -----------------------------
# Load coverage
# -----------------------------
coverage = pd.read_csv(primary_cov, sep="\t")
coverage = coverage.rename(columns={
    "#rname": "reference",
    "coverage": "breadth_percent",
    "meandepth": "mean_depth"
})

coverage["breadth"] = coverage["breadth_percent"] / 100

# -----------------------------
# Merge all
# -----------------------------
df = primary.merge(strict, on="reference")
df = df.merge(coverage[["reference", "breadth", "mean_depth"]], on="reference")

# -----------------------------
# Normalize by length (per kb)
# -----------------------------
df["primary_per_kb"] = df["primary_reads"] / (df["length"] / 1000)
df["unique_per_kb"] = df["unique_reads"] / (df["length"] / 1000)

# -----------------------------
# Scoring formula
# -----------------------------
df["score_base"] = (2 * df["unique_per_kb"]) + df["primary_per_kb"]
df["score"] = df["score_base"] * (2* df["breadth"])

# -----------------------------
# Rank references
# -----------------------------
df = df.sort_values("score", ascending=False).reset_index(drop=True)

# Confidence ratio
if len(df) > 1 and df.loc[1, "score"] > 0:
    ratio = df.loc[0, "score"] / df.loc[1, "score"]
else:
    ratio = float("inf")

if ratio >= 2:
    confidence = "strong"
elif ratio >= 1.3:
    confidence = "moderate"
else:
    confidence = "ambiguous"

df["confidence"] = ""
df.loc[0, "confidence"] = confidence

# Save full ranking
df.to_csv(output_file, sep="\t", index=False)