# Cromwell Setup - Offline VRS Workflows

## Quick Start

**Prerequisites:**
- Java (JRE/JDK 11+)
- Docker daemon running
- Pre-built `cerfac:vrs-offline` Docker image

**Run workflow:**
```bash
java -Dconfig.file=workflows/combined_gnomad_clinvar/cromwell.conf \
  -jar tools/cromwell-89.jar run \
  workflows/combined_gnomad_clinvar/merge_clinical_functional_data_vrs_offline.wdl \
  --inputs workflows/test/merge_clinical_functional_data_vrs_offline.brca1.input.json
```

**Verify success:**
```bash
# Look for "workflow finished with status 'Succeeded'"
# Check outputs:
ls cromwell-executions/merge_clinical_data/*/call-merge_vrs_files/execution/BRCA1_variants_functional_clinical.csv
```

---

## Installation

**Cromwell**: Pre-downloaded in `tools/cromwell-89.jar` (250 MB)

**Java:**
```bash
# macOS
brew install java

# Linux (Ubuntu/Debian)
sudo apt-get install default-jre
```

**Docker:**
```bash
# Ensure Docker daemon is running
docker --version
```

## Setup

### 1. Build Offline VRS Docker Image

```bash
cd workflows/compute_vrs_digests
docker build -t cerfac:vrs-offline .
```
- **Time**: 10-15 minutes
- **Size**: 13.3 GB (SeqRepo GRCh38 + UTA + PostgreSQL)
- One-time setup

### 2. Configuration

`workflows/combined_gnomad_clinvar/cromwell.conf` is pre-configured:
- Backend: Local (runs on this machine)
- Docker: Enabled with proper volume mounts
- VRS tasks: Automatically use `cerfac:vrs-offline` image

No changes needed for basic use.

---

## Running Workflows

### Offline VRS Workflow

Computes VRS digests for three data sources and merges results on VRS ID:

```bash
java -Dconfig.file=workflows/combined_gnomad_clinvar/cromwell.conf \
  -jar tools/cromwell-89.jar run \
  workflows/combined_gnomad_clinvar/merge_clinical_functional_data_vrs_offline.wdl \
  --inputs workflows/test/merge_clinical_functional_data_vrs_offline.brca1.input.json
```

**Input JSON format:**
```json
{
  "merge_clinical_data.GENE_NAME": "BRCA1",
  "merge_clinical_data.VARIANTS_FILE": "/path/to/gnomad_variants.csv",
  "merge_clinical_data.FUNCTIONAL_SCORES": "/path/to/functional_scores.csv",
  "merge_clinical_data.CLINICAL_DATA": "/path/to/clinical_data.csv"
}
```

**Output:**
- Location: `cromwell-executions/merge_clinical_data/[workflow-id]/call-merge_vrs_files/execution/`
- File: `[GENE_NAME]_variants_functional_clinical.csv`
- Contains: gnomAD + functional scores + clinical data, merged on VRS IDs

**Performance:**
- Per-variant: 0.3-0.5 seconds
- 3,784 variants (BRCA1 test): ~30 minutes
- Completely offline (zero network calls)

### Validate Syntax

```bash
java -jar tools/cromwell-89.jar validate \
  workflows/combined_gnomad_clinvar/merge_clinical_functional_data_vrs_offline.wdl
```

---

## What Gets Computed

The workflow runs three parallel VRS tasks, then merges:

| Task | Input | Output | Format |
|------|-------|--------|--------|
| vrs_variants | gnomAD CSV | VRS digests + coordinates | chr-pos-ref-alt |
| vrs_scores | Functional assay CSV | VRS digests + amino acid info | chr/pos/ref/alt columns |
| vrs_clinical | Clinical data CSV | VRS digests | Variant IDs |
| merge_vrs_files | All three above | Unified dataset | All columns merged on vrs_id |

**Output columns include:**
- gnomAD: population frequencies, variant annotations
- Functional scores: assay scores, amino acid changes
- Clinical: BRIDGES, CARRIERS, UKB odds ratios and case/control counts
- VRS: vrs_id, vrs_digest, genomic_hgvs, vrs_start, vrs_end

---

## Offline System Details

### Docker Image Contents

- **SeqRepo GRCh38**: Local sequence reference (no remote queries)
- **UTA Database**: Transcript annotations (embedded PostgreSQL)
- **ga4gh-vrs**: VRS computation library
- **Environment**: All variables pre-configured for offline operation

### Environment Variables (Auto-configured in image)

```bash
GA4GH_VRS_DATAPROXY_URI=seqrepo+file:///seqrepo-GRCh38/master
UTA_DB_URL=postgresql://postgres:postgres@localhost/uta/uta_20241220
BIOUTILS_NCBI_RETRIES=0
```

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| "Java not found" | Install JRE: `brew install java` |
| "Docker daemon not running" | Start Docker application |
| "Cannot find cerfac:vrs-offline" | Build image: `docker build -t cerfac:vrs-offline workflows/compute_vrs_digests/` |
| Task timeout | Increase memory in cromwell.conf: `--memory=16g` |
| Out of memory (OOM) | Increase memory: `--memory=16g` for larger variant files |

## Output Locations

```
cromwell-executions/merge_clinical_data/
├── [workflow-id]/
│   ├── call-vrs_variants/execution/variants_with_vrs.csv
│   ├── call-vrs_scores/execution/functional_scores_with_vrs.csv
│   ├── call-vrs_clinical/execution/clinical_data_with_vrs.csv
│   └── call-merge_vrs_files/execution/[GENE]_variants_functional_clinical.csv
```

View workflow logs:
```bash
tail -100 cromwell-executions/merge_clinical_data/[workflow-id]/*/execution/stdout
```

---

## References

- [Cromwell](https://cromwell.readthedocs.io/)
- [WDL](https://github.com/openwdl/wdl)
- [GA4GH VRS](https://vrs.ga4gh.org/)

