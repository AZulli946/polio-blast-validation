# Poliovirus BLAST Validation Pipeline

Extracts poliovirus-classified reads from EsViritu metagenomics BAM files, validates them via BLAST against a curated database of complete Enterovirus C genomes that explicitly includes all poliovirus genomes, and annotates results with genome region mapping.

## Quick Start

```bash
cd polio-blast-validation
pip install -r requirements.txt    # pandas, pysam, google-cloud-storage, biopython, pyyaml
python run_pipeline.py             # runs all 5 steps
```

## Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 1 | `identify_hits.py` | Find poliovirus hits in precomputed dashboard data |
| 2 | `extract_reads.py` | Download BAMs from GCS, extract poliovirus-aligned reads |
| 3 | `build_blastdb.py` | Download complete Enterovirus C genomes from NCBI, build BLAST DB |
| 4 | `blast_validate.py` | Run blastn, join results with read metadata |
| 5 | `annotate_regions.py` | Map BLAST coordinates to genome regions (VP1, 3D, etc.) |

Run individual steps: `python run_pipeline.py --step identify`

## Requirements

- Python 3.9+
- BLAST+ (`blastn`, `makeblastdb`) on PATH
- GCS service account key at `../shiny_dashboard/gcs_service_account.json`
- Precomputed virus data at `../shiny_dashboard/dashboard_data/precomputed_virus.tsv.gz`

## Output

Final table: `results/blast_validated.tsv` with columns:
- `read_id`, `site`, `collection_date`, `esviritu_classification`
- `top_blast_hit`, `percent_identity`, `query_length`, `alignment_length`
- `subject_start`, `subject_end`, `genome_region`
- `evalue`, `bitscore`, `sample_ID`, `delivery_date`, `esviritu_accession`

Reads with no BLAST hit are included with empty BLAST columns.

## Configuration

Edit `config.yaml` to change paths, subspecies patterns, BLAST parameters, or the Enterovirus C Entrez query used to build the reference database. Step 5 is selected by `annotate_module` (here `annotate_regions`); the build-DB scope summary is driven by `reference_classes`.

## Shared core with newcastle-blast-validation

This repo shares an identical "core" (steps 1–4 + the diagram) with the sibling
`newcastle-blast-validation` repo. Only the following files are target-specific:
`config.yaml`, `README.md`, `genome_regions.tsv`, and the step-5 module
`scripts/annotate_regions.py`.

The shared files must stay byte-identical across both repos. After editing any of
them, verify with:

```bash
python tools/check_sync.py --other ../newcastle-blast-validation
```

Shared files: `run_pipeline.py`, `requirements.txt`, `.gitignore`,
`scripts/{identify_hits,extract_reads,build_blastdb,blast_validate,make_pipeline_diagram}.py`,
`tools/check_sync.py`.
