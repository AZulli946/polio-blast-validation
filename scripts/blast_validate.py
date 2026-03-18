"""Step 4: Run BLAST validation and join results with read metadata."""

import os
import shutil
import subprocess
import tempfile
import pandas as pd


BLAST_COLUMNS = [
    "qseqid", "sseqid", "stitle", "pident", "length",
    "qlen", "sstart", "send", "slen", "evalue", "bitscore",
]


def run(config, data_dir, results_dir):
    db_name = config["blast_db_name"]
    query_path = results_dir / "poliovirus_reads.fasta"
    blast_out = results_dir / "blast_results.tsv"
    meta_path = results_dir / "read_metadata.tsv"

    if not query_path.exists() or query_path.stat().st_size == 0:
        print("No reads to BLAST (empty FASTA). Skipping.")
        pd.DataFrame().to_csv(results_dir / "blast_validated.tsv", sep="\t", index=False)
        return

    db_files = list(data_dir.glob(f"blastdb/{db_name}.*"))
    if not db_files:
        raise FileNotFoundError(f"BLAST database not found. Run build_blastdb first.")

    max_target_seqs = config.get("max_target_seqs", 1)
    evalue = config.get("evalue", 1e-5)

    outfmt = "6 " + " ".join(BLAST_COLUMNS)

    # BLAST+ can't handle paths with spaces, so copy files to a temp dir
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_query = os.path.join(tmpdir, "query.fasta")
        tmp_out = os.path.join(tmpdir, "blast_results.tsv")

        shutil.copy2(str(query_path), tmp_query)

        # Copy BLAST DB files to tmpdir
        db_dir = data_dir / "blastdb"
        for f in db_dir.iterdir():
            if f.name.startswith(db_name):
                shutil.copy2(str(f), os.path.join(tmpdir, f.name))

        tmp_db = os.path.join(tmpdir, db_name)

        cmd = [
            "blastn",
            "-query", tmp_query,
            "-db", tmp_db,
            "-outfmt", outfmt,
            "-max_target_seqs", str(max_target_seqs),
            "-evalue", str(evalue),
            "-out", tmp_out,
        ]

        print(f"Running BLAST ({max_target_seqs} max target seqs, evalue {evalue})...")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"STDERR: {result.stderr}")
            raise RuntimeError(f"blastn failed with exit code {result.returncode}")

        # Copy results back
        shutil.copy2(tmp_out, str(blast_out))

    # Parse BLAST results
    if blast_out.stat().st_size == 0:
        print("WARNING: BLAST produced no hits.")
        blast_df = pd.DataFrame(columns=BLAST_COLUMNS)
    else:
        blast_df = pd.read_csv(blast_out, sep="\t", header=None, names=BLAST_COLUMNS)

    print(f"  BLAST hits: {len(blast_df):,}")
    print(f"  Unique queries with hits: {blast_df['qseqid'].nunique()}")

    # Load read metadata
    meta_df = pd.read_csv(meta_path, sep="\t", dtype=str)

    # Join BLAST results to metadata on unique read ID (left join to keep reads with no hit)
    merged = meta_df.merge(blast_df, left_on="read_uid", right_on="qseqid", how="left")

    if "qseqid" in merged.columns:
        merged.drop(columns=["qseqid"], inplace=True)

    n_no_hit = merged["sseqid"].isna().sum()
    n_with_hit = merged["sseqid"].notna().sum()
    print(f"  Reads with BLAST hit: {n_with_hit}")
    print(f"  Reads without BLAST hit: {n_no_hit}")

    merged.to_csv(results_dir / "blast_merged.tsv", sep="\t", index=False)
    print(f"  Merged results written to {results_dir / 'blast_merged.tsv'}")

    return merged
