"""Step 3: Download complete target genomes from NCBI and build a BLAST database.

Shared core (identical across repos). Two config-driven extension points:
  - ``reference_classes``: list of {label, patterns} used to summarize DB scope.
  - ``reference_annotator``: optional module name (in scripts/) exposing
    ``annotation_fields`` + ``find_annotation(record, config)``. When set, a
    per-reference annotation TSV is written for downstream steps to consume.
"""

import csv
import importlib
import io
import os
import shutil
import subprocess
import tempfile
from collections import Counter

from Bio import Entrez, SeqIO


def _classify_reference(description, reference_classes):
    """Return the first class label whose patterns match the description."""
    desc = description.lower()
    for cls in reference_classes:
        for pattern in cls.get("patterns", []):
            if pattern.lower() in desc:
                return cls["label"]
    return "other"


def _load_annotator(config):
    """Return the configured reference-annotator module, or None."""
    module_name = config.get("reference_annotator")
    if not module_name:
        return None
    return importlib.import_module(f"scripts.{module_name}")


def run(config, data_dir):
    db_name = config["blast_db_name"]
    entrez_query = config["entrez_query"]
    entrez_email = config["entrez_email"]
    target_name = config.get("target_name", "target virus")
    batch_size = int(config.get("entrez_batch_size", 200))
    reference_classes = config.get("reference_classes", [])
    annotator = _load_annotator(config)

    Entrez.email = entrez_email

    db_dir = data_dir / "blastdb"
    db_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = db_dir / "genomes.fasta"
    genbank_path = db_dir / "references.gb"
    annotation_path = db_dir / "reference_annotations.tsv"

    # Search NCBI for complete target-virus genomes.
    print(f"Searching NCBI: {entrez_query}")
    handle = Entrez.esearch(db="nucleotide", term=entrez_query, retmax=10000)
    record = Entrez.read(handle)
    handle.close()

    ids = record["IdList"]
    count = int(record["Count"])
    print(f"  Found {count} sequences, fetching {len(ids)}")

    if not ids:
        raise RuntimeError("No sequences found. Check your Entrez query.")

    print(f"Downloading GenBank records for {target_name}...")
    counts = Counter()
    annotation_rows = []
    n_seqs = 0

    with open(genbank_path, "w") as genbank_out, open(fasta_path, "w") as fasta_out:
        for start in range(0, len(ids), batch_size):
            batch_ids = ids[start:start + batch_size]
            print(f"  Fetching {start + 1}-{start + len(batch_ids)} of {len(ids)}")
            handle = Entrez.efetch(
                db="nucleotide",
                id=",".join(batch_ids),
                rettype="gbwithparts",
                retmode="text",
            )
            genbank_data = handle.read()
            handle.close()
            genbank_out.write(genbank_data)

            for seq_record in SeqIO.parse(io.StringIO(genbank_data), "genbank"):
                n_seqs += 1
                counts[_classify_reference(seq_record.description, reference_classes)] += 1
                SeqIO.write(seq_record, fasta_out, "fasta")

                if annotator:
                    annotation = annotator.find_annotation(seq_record, config)
                    if annotation:
                        annotation_rows.append(annotation)

    print(f"  Written to {fasta_path}")
    print(f"  GenBank records written to {genbank_path}")

    if annotator:
        with open(annotation_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=annotator.annotation_fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(annotation_rows)
        print(f"  Reference annotations written to {annotation_path}")

    # Summarize downloaded references so the DB scope is obvious in logs.
    print(f"  Total sequences in FASTA: {n_seqs}")
    for cls in reference_classes:
        label = cls["label"]
        if counts[label]:
            print(f"    {label}: {counts[label]}")
    if counts["other"]:
        print(f"    other: {counts['other']}")
    if annotator:
        print(f"  References with usable annotation: {len(annotation_rows)}")
        if not annotation_rows:
            print("  WARNING: No usable annotations found; downstream annotation may be unresolved.")

    # Build BLAST database in a temp dir (BLAST+ can't handle paths with spaces)
    print("Building BLAST database...")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_fasta = os.path.join(tmpdir, "genomes.fasta")
        tmp_db = os.path.join(tmpdir, db_name)
        shutil.copy2(str(fasta_path), tmp_fasta)

        cmd = [
            "makeblastdb",
            "-in", tmp_fasta,
            "-dbtype", "nucl",
            "-out", tmp_db,
            "-parse_seqids",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"STDERR: {result.stderr}")
            raise RuntimeError(f"makeblastdb failed with exit code {result.returncode}")

        print(result.stdout)

        # Copy DB files back to the actual data directory
        for f in os.listdir(tmpdir):
            if f.startswith(db_name) and f != "genomes.fasta":
                shutil.copy2(os.path.join(tmpdir, f), str(db_dir / f))

    print(f"  BLAST database built at {db_dir / db_name}")
    return db_dir / db_name
