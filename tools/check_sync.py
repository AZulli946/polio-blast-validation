#!/usr/bin/env python3
"""Verify the shared "core" pipeline files are byte-identical to a sibling repo.

The newcastle and polio BLAST-validation repos deliberately share an identical
set of step 1-4 + diagram files. This tool diffs them so drift is caught early.

Usage:
    python tools/check_sync.py [--other PATH]

--other defaults to auto-detecting a sibling "*-blast-validation" directory next
to this repo. Exits non-zero if any shared file differs or is missing.
"""

import argparse
import sys
from pathlib import Path

# Files that MUST be byte-identical across the newcastle and polio repos.
# Target-specific files are intentionally excluded:
#   config.yaml, README.md, genome_regions.tsv,
#   scripts/{annotate_regions,classify_fusion_junction,fusion_rules,fusion_annotator}.py
SHARED_FILES = [
    "run_pipeline.py",
    "requirements.txt",
    ".gitignore",
    "scripts/identify_hits.py",
    "scripts/extract_reads.py",
    "scripts/build_blastdb.py",
    "scripts/blast_validate.py",
    "scripts/make_pipeline_diagram.py",
    "tools/check_sync.py",
]


def find_sibling(repo_dir):
    candidates = [
        p for p in repo_dir.parent.iterdir()
        if p.is_dir() and p.name.endswith("-blast-validation") and p != repo_dir
    ]
    return candidates[0] if len(candidates) == 1 else None


def main():
    repo_dir = Path(__file__).resolve().parent.parent

    parser = argparse.ArgumentParser(description="Check shared-core files are in sync")
    parser.add_argument("--other", help="Path to the sibling repo to compare against")
    args = parser.parse_args()

    other = Path(args.other).resolve() if args.other else find_sibling(repo_dir)
    if other is None or not other.is_dir():
        print("ERROR: could not locate sibling repo. Pass --other PATH.")
        return 2

    print(f"Comparing:\n  this:  {repo_dir}\n  other: {other}\n")

    n_diff = 0
    for rel in SHARED_FILES:
        a, b = repo_dir / rel, other / rel
        if not a.exists():
            print(f"MISSING (this)   {rel}")
            n_diff += 1
        elif not b.exists():
            print(f"MISSING (other)  {rel}")
            n_diff += 1
        elif a.read_bytes() != b.read_bytes():
            print(f"DIFFER           {rel}")
            n_diff += 1
        else:
            print(f"IDENTICAL        {rel}")

    print()
    if n_diff:
        print(f"OUT OF SYNC: {n_diff} file(s) differ. Reconcile before committing.")
        return 1
    print("All shared-core files are in sync.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
