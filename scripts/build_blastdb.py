"""Step 3: Download complete Enterovirus C genomes and build a BLAST database."""

import os
import shutil
import subprocess
import tempfile
from collections import Counter
from Bio import Entrez, SeqIO


def _classify_reference(description):
    desc = description.lower()
    if "poliovirus" in desc:
        return "poliovirus"
    if "coxsackievirus" in desc:
        return "coxsackievirus"
    if "enterovirus c" in desc or "human enterovirus" in desc:
        return "enterovirus_c_other"
    return "other"


def run(config, data_dir):
    db_name = config["blast_db_name"]
    entrez_query = config["entrez_query"]
    entrez_email = config["entrez_email"]

    Entrez.email = entrez_email

    db_dir = data_dir / "blastdb"
    db_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = db_dir / "genomes.fasta"

    # Search NCBI for complete Enterovirus C genomes, explicitly retaining poliovirus.
    print(f"Searching NCBI: {entrez_query}")
    handle = Entrez.esearch(db="nucleotide", term=entrez_query, retmax=10000)
    record = Entrez.read(handle)
    handle.close()

    ids = record["IdList"]
    count = int(record["Count"])
    print(f"  Found {count} sequences, fetching {len(ids)}")

    if not ids:
        raise RuntimeError("No sequences found. Check your Entrez query.")

    # Fetch sequences in FASTA format
    print("Downloading sequences...")
    handle = Entrez.efetch(db="nucleotide", id=ids, rettype="fasta", retmode="text")
    fasta_data = handle.read()
    handle.close()

    fasta_path.write_text(fasta_data)
    print(f"  Written to {fasta_path}")

    # Summarize downloaded references so the DB scope is obvious in logs.
    counts = Counter()
    n_seqs = 0
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        n_seqs += 1
        counts[_classify_reference(record.description)] += 1
    print(f"  Total sequences in FASTA: {n_seqs}")
    for label in ("poliovirus", "coxsackievirus", "enterovirus_c_other", "other"):
        if counts[label]:
            print(f"    {label}: {counts[label]}")

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
