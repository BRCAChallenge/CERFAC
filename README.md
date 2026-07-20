# CERFAC

Cloud-based, Evidence-based Rigorous Functional Assay Calibration

A bioinformatics pipeline for validating functional assay predictions against clinical observational data using WDL workflows on Terra or locally.

## Quick Start

1. **New to CERFAC?** Read [CERFAC_documentation.md](CERFAC_documentation.md) for the full step-by-step guide
2. **Setting up locally?** See [CROMWELL_SETUP.md](CROMWELL_SETUP.md) for workflow execution
3. **Want to understand the architecture?** Check [CLAUDE.md](CLAUDE.md) for detailed technical information

## What It Does

CERFAC combines three data sources to calibrate functional assay predictions:

- **ClinVar**: Clinical variant classification from NCBI
- **gnomAD**: Population variant frequency data
- **Your assay scores**: Functional predictions you provide

The pipeline normalizes variants across formats, merges the data, and performs statistical calibration (OR calculations) via Jupyter notebook.

## Running on Terra (Production)

See [CERFAC_documentation.md](CERFAC_documentation.md) — complete step-by-step guide for cloning a workspace and running workflows.

## Running Locally (Development)

See [CROMWELL_SETUP.md](CROMWELL_SETUP.md) for workflow execution or [CLAUDE.md](CLAUDE.md) for technical details.

## Architecture

```
ClinVar Query → ┐
                ├→ Merge & Normalize → Statistical Analysis
gnomAD Query  → ┤   (VRS API)          (R Jupyter notebook)
Your Data ────→ ┘
```

Three Docker images handle different stages:
- `cerfac-clinvar` — NCBI EDirect extraction
- `cerfac-merge` — Data merging and API calls
- `cerfac-gnomad` — Population frequency lookup via Hail

## Documentation

| Document | Purpose |
|----------|---------|
| [CERFAC_documentation.md](CERFAC_documentation.md) | User guide for Terra |
| [CLAUDE.md](CLAUDE.md) | Developer reference |
| [CROMWELL_SETUP.md](CROMWELL_SETUP.md) | Local workflow execution |
| [docker/README.md](docker/README.md) | Docker image building |

## Key Features

- **Variant normalization**: Handles HGVS coding, HGVS genomic, and gnomAD VCF formats
- **Scalable**: Runs on Terra cloud or local machine via Cromwell
- **Reproducible**: WDL workflows with containerized dependencies
- **Flexible**: Accepts TSV, CSV, or TXT input files

## Requirements

- Docker (for containerized task execution)
- Java Runtime (for Cromwell workflow engine)
- Sufficient disk space (~2 GB for Docker images)

## Questions?

Refer to the documentation files above or check [CLAUDE.md](CLAUDE.md) for technical details.
