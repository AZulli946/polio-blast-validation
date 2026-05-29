"""Step 1: Identify target-virus hits in precomputed virus data."""

import pandas as pd


def run(config, results_dir):
    virus_path = config["precomputed_virus_data"]
    patterns = config["subspecies_patterns"]
    target_name = config.get("target_name", "target virus")
    target_slug = config.get("target_slug", "target")

    print(f"Reading precomputed virus data from {virus_path}")
    df = pd.read_csv(virus_path, sep="\t", compression="gzip", dtype=str)
    print(f"  Total rows: {len(df):,}")

    # Build combined regex from all patterns (non-capturing groups to avoid warning)
    combined = "|".join(f"(?:{p})" for p in patterns)
    mask = df["subspecies"].str.contains(combined, regex=True, na=False)
    hits = df.loc[mask].copy()
    print(f"  {target_name} hits: {len(hits):,}")

    delivery_filter = config.get("delivery_filter") or []
    if delivery_filter:
        delivery_filter = [str(d) for d in delivery_filter]
        hits = hits.loc[hits["delivery_date"].astype(str).isin(delivery_filter)].copy()
        print(f"  After delivery filter {delivery_filter}: {len(hits):,}")

    if hits.empty:
        print(f"WARNING: No {target_name} hits found.")

    # Select and output relevant columns
    keep_cols = [
        "sample_ID", "Accession", "Name", "subspecies", "City", "Date",
        "delivery_date", "read_count", "reads_per_million",
    ]
    # Only keep columns that exist
    keep_cols = [c for c in keep_cols if c in hits.columns]
    hits = hits[keep_cols]

    out_path = results_dir / f"{target_slug}_hits.tsv"
    hits.to_csv(out_path, sep="\t", index=False)
    print(f"  Written to {out_path}")
    print(f"  Unique samples: {hits['sample_ID'].nunique()}")
    print(f"  Unique subspecies: {hits['subspecies'].unique().tolist()}")

    return hits
