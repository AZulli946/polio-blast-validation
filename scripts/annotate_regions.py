"""Step 5: Map BLAST alignment coordinates to poliovirus genome regions."""

import re
import pandas as pd


def _load_genome_regions(regions_path):
    """Load genome_regions.tsv into a dict keyed by serotype."""
    df = pd.read_csv(regions_path, sep="\t")
    regions = {}
    for serotype, group in df.groupby("serotype"):
        regions[serotype] = list(group[["region", "start", "end"]].itertuples(index=False, name=None))
    return regions


def _determine_serotype(stitle):
    """Parse serotype from BLAST subject title."""
    if pd.isna(stitle):
        return None

    stitle_lower = stitle.lower()

    # Match patterns like "poliovirus 1", "poliovirus type 1", "Human poliovirus 1"
    m = re.search(r'poliovirus\s*(?:type\s*)?(\d)', stitle_lower)
    if m:
        return f"PV{m.group(1)}"

    # Match "PV1", "PV2", "PV3"
    m = re.search(r'\bpv(\d)\b', stitle_lower)
    if m:
        return f"PV{m.group(1)}"

    # Match "Sabin 1/2/3" (vaccine strains)
    m = re.search(r'sabin[- ]*(\d)', stitle_lower)
    if m:
        return f"PV{m.group(1)}"

    return None


def _map_to_regions(sstart, send, region_list):
    """Map alignment coordinates to genome regions.

    Returns comma-separated region names that overlap with [sstart, send].
    """
    if pd.isna(sstart) or pd.isna(send):
        return ""

    sstart = int(float(sstart))
    send = int(float(send))

    # Handle reverse complement alignments
    if sstart > send:
        sstart, send = send, sstart

    overlapping = []
    for region_name, region_start, region_end in region_list:
        # Check overlap
        if sstart <= region_end and send >= region_start:
            overlapping.append(region_name)

    return ", ".join(overlapping) if overlapping else "intergenic/unknown"


def run(config, project_dir, results_dir):
    regions_path = project_dir / "genome_regions.tsv"
    merged_path = results_dir / "blast_merged.tsv"

    if not merged_path.exists():
        print("No blast_merged.tsv found. Skipping annotation.")
        return

    regions = _load_genome_regions(regions_path)
    print(f"Loaded genome regions for serotypes: {list(regions.keys())}")

    df = pd.read_csv(merged_path, sep="\t", dtype=str)

    if df.empty:
        print("No data to annotate.")
        df.to_csv(results_dir / "blast_validated.tsv", sep="\t", index=False)
        return

    # Determine serotype for each hit
    df["serotype"] = df["stitle"].apply(_determine_serotype)

    # Default to PV1 if serotype can't be determined
    default_serotype = "PV1"
    n_default = df["serotype"].isna().sum() - df["stitle"].isna().sum()  # exclude no-hit rows
    if n_default > 0:
        print(f"  WARNING: {n_default} hits with undetermined serotype, defaulting to {default_serotype}")
    df["serotype"] = df["serotype"].fillna(default_serotype)

    # Map to genome regions
    def annotate_row(row):
        serotype = row["serotype"]
        region_list = regions.get(serotype, regions.get(default_serotype, []))
        return _map_to_regions(row.get("sstart"), row.get("send"), region_list)

    df["genome_region"] = df.apply(annotate_row, axis=1)

    # Drop metadata query_length (redundant with BLAST qlen)
    if "query_length" in df.columns and "qlen" in df.columns:
        df.drop(columns=["query_length"], inplace=True)

    # Build final output table with specified columns
    output_cols = {
        "read_id": "read_id",
        "City": "site",
        "Date": "collection_date",
        "subspecies": "esviritu_classification",
        "stitle": "top_blast_hit",
        "pident": "percent_identity",
        "qlen": "query_length",
        "length": "alignment_length",
        "sstart": "subject_start",
        "send": "subject_end",
        "genome_region": "genome_region",
        "evalue": "evalue",
        "bitscore": "bitscore",
        "sample_ID": "sample_ID",
        "delivery_date": "delivery_date",
        "accession": "esviritu_accession",
        "read_sequence": "read_sequence",
    }

    # Only rename columns that exist
    rename_map = {k: v for k, v in output_cols.items() if k in df.columns}
    final = df.rename(columns=rename_map)

    # Select output columns (in order, if they exist)
    desired_order = list(output_cols.values())
    final_cols = [c for c in desired_order if c in final.columns]
    final = final[final_cols]

    out_path = results_dir / "blast_validated.tsv"
    final.to_csv(out_path, sep="\t", index=False)

    # Summary stats
    n_total = len(final)
    n_validated = final["top_blast_hit"].notna().sum()
    n_unvalidated = final["top_blast_hit"].isna().sum()

    print(f"\nFinal results: {n_total} reads")
    print(f"  Validated (BLAST hit): {n_validated}")
    print(f"  Unvalidated (no BLAST hit): {n_unvalidated}")

    if n_validated > 0:
        validated = final[final["top_blast_hit"].notna()]
        print(f"\n  Genome region distribution:")
        region_counts = validated["genome_region"].value_counts()
        for region, count in region_counts.items():
            print(f"    {region}: {count}")

        print(f"\n  Top BLAST hits:")
        hit_counts = validated["top_blast_hit"].value_counts().head(10)
        for hit, count in hit_counts.items():
            # Truncate long titles
            display = hit[:80] + "..." if len(str(hit)) > 80 else hit
            print(f"    {display}: {count}")

    print(f"\nWritten to {out_path}")

    return final
