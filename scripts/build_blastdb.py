"""Step 3: Download complete poliovirus genomes and build BLAST database."""

import os
import shutil
import subprocess
import tempfile
from Bio import Entrez, SeqIO


def run(config, data_dir):
    db_name = config["blast_db_name"]
    entrez_query = config["entrez_query"]
    entrez_email = config["entrez_email"]

    Entrez.email = entrez_email

    db_dir = data_dir / "blastdb"
    db_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = db_dir / "genomes.fasta"

    # Search NCBI for poliovirus complete genomes
    print(f"Searching NCBI: {entrez_query}")
    handle = Entrez.esearch(db="nucleotide", term=entrez_query, retmax=500)
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

    # Count sequences
    n_seqs = sum(1 for _ in SeqIO.parse(str(fasta_path), "fasta"))
    print(f"  Total sequences in FASTA: {n_seqs}")

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
