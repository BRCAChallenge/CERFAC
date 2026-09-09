# VRS Digest Computation & Merge Workflow

Computes GA4GH Variation Representation Specification (VRS) digests for genomic variants and merges multiple variant datasets using VRS IDs as the common key.

## Overview

**Three-stage workflow:**

1. **VRS Digest Computation**: Converts variant data from multiple formats (VCF, genomic coords, etc.) into standardized VRS digests
2. **Data Integration**: Merges gnomAD reference variants with user-provided clinical and functional assay data
3. **Unified Output**: Single CSV with all data for a gene, keyed by VRS digest

**Input data types supported:**
- **gnomAD variants**: Pre-computed VCF format variants (via get_gnomad_variants workflow)
- **Functional assay scores**: Separate chr/pos/ref/alt columns
- **Clinical case-control data**: Variant IDs in VCF format

**Output:**
- Three individual CSVs with VRS digests (gnomAD, functional, clinical)
- One merged CSV joining all three datasets by VRS ID

## Workflow Architecture

```
User Input:
  - Gene name (e.g., BRCA1)
  - gnomAD variants file (required)
  - Functional assay file (optional)
  - Clinical data file (optional)
         |
         v
   [Compute VRS Digests]
    - gnomAD → VRS IDs (chr17-43045682-T-C → ga4gh:VA....)
    - Functional → VRS IDs (chr/pos/ref/alt → ga4gh:VA....)
    - Clinical → VRS IDs (13-32314431-T-C → ga4gh:VA....)
         |
         v
   [Merge by VRS ID]
    - Outer join on vrs_id
    - Retain all variants from all sources
    - Columns prefixed by source (gnomad_*, functional_*, clinical_*)
         |
         v
   Output:
    - BRCA1_merged_vrs.csv (unified dataset)
```

## Files

- `compute_vrs_digests.py` — Batch variant processor using AlleleTranslator
- `compute_vrs_digests.wdl` — Cromwell workflow (4 tasks: 3 VRS + 1 merge)
- `Dockerfile` — Extends eugene75 vrs-in-a-box with our script
- `README.md` — This file
- `test_inputs.brca1.json` — Example with all three data types

## Usage

### Docker (Local Testing)

**Single data type:**
```bash
docker run --rm \
  -v /path/to/data:/data \
  cerfac/vrs-digests:grch38 \
  python3 /vrs-python/compute_vrs_digests.py \
    --input /data/variants.csv \
    --output /data/variants_with_vrs.csv \
    --variant-type gnomad
```

Variant types: `gnomad`, `functional_assay`, `clinical`

### Cromwell Workflow (Includes Merge)

```bash
java -jar cromwell.jar run compute_vrs_digests.wdl \
  --inputs test_inputs.brca1.json
```

Outputs:
- `BRCA1_gnomad_with_vrs.csv` (gnomAD + VRS digests)
- `BRCA1_functional_assay_with_vrs.csv` (Functional + VRS digests)
- `BRCA1_clinical_with_vrs.csv` (Clinical + VRS digests)
- `BRCA1_merged_vrs.csv` (All three merged by VRS ID)

## VRS Digest Computation

Converts variants to standardized VRS format using ga4gh-vrs AlleleTranslator:

**Input → VRS ID mapping:**
- gnomAD VCF (`chr17-43045682-T-C`) → `ga4gh:VA.9n9ADvjO6BOKLPQvP7lTVScJPhghON5q`
- Functional coords (`17:43045682:T:C`) → `ga4gh:VA.9n9ADvjO6BOKLPQvP7lTVScJPhghON5q`
- Clinical VCF (`17-43045682-T-C`) → `ga4gh:VA.9n9ADvjO6BOKLPQvP7lTVScJPhghON5q`

**Same variant = same digest** → merge key works across all data sources

## Merge Strategy

**Join type:** Outer join (retains all variants from all sources)

**Merge key:** `vrs_id` (VRS digest common to all data types)

**Column naming:** Prefixed by source to avoid conflicts
- `gnomad_vrs_digest`, `gnomad_genomic_hgvs`
- `functional_vrs_digest`, `functional_genomic_hgvs`
- `clinical_vrs_digest`, `clinical_genomic_hgvs`
- Plus all original columns from each source

**Result:** One row per unique variant, with NULL values where a variant appears in some but not all data sources.

## Building the Docker Image

```bash
docker build -t brcachallenge/cerfac:vrs-digests-latest .
```

Uses SeqRepo and vrs-python from eugene75 base image. Published to `brcachallenge/cerfac:vrs-digests-latest` on Docker Hub.

## Testing

Tested with 50-variant subsets:
- gnomAD variants: **100% success** (50/50)
- Functional assay: **100% success** (50/50)
- Clinical data: **100% success** (50/50)
- Merge operation: Verified VRS ID matching across sources

Processing time: ~2-3 minutes per 100 variants

## References

- [GA4GH VRS Specification](https://vrs.ga4gh.org/)
- [vrs-python](https://github.com/ga4gh/vrs-python)
- [vrs-in-a-box](https://github.com/ehclark/vrs-in-a-box)
