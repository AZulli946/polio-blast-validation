"""Step 1: Identify poliovirus hits in precomputed virus data."""

import re
import pandas as pd


def run(config, results_dir):
    virus_path = config["precomputed_virus_data"]
    patterns = config["subspecies_patterns"]

    print(f"Reading precomputed virus data from {virus_path}")
    df = pd.read_csv(virus_path, sep="\t", compression="gzip", dtype=str)
    print(f"  Total rows: {len(df):,}")

    # Build combined regex from all patterns (non-capturing groups to avoid warning)
    combined = "|".join(f"(?:{p})" for p in patterns)
    mask = df["subspecies"].str.contains(combined, regex=True, na=False)
    hits = df.loc[mask].copy()
    print(f"  Poliovirus hits: {len(hits):,}")

    if hits.empty:
        print("WARNING: No poliovirus hits found.")

    # Select and output relevant columns
    keep_cols = [
        "sample_ID", "Accession", "Name", "subspecies", "City", "Date",
        "delivery_date", "read_count", "reads_per_million",
    ]
    # Only keep columns that exist
    keep_cols = [c for c in keep_cols if c in hits.columns]
    hits = hits[keep_cols]

    out_path = results_dir / "poliovirus_hits.tsv"
    hits.to_csv(out_path, sep="\t", index=False)
    print(f"  Written to {out_path}")
    print(f"  Unique samples: {hits['sample_ID'].nunique()}")
    print(f"  Unique subspecies: {hits['subspecies'].unique().tolist()}")

    return hits
