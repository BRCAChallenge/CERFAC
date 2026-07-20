# Example Data Files for CERFAC

This directory contains example data files for testing and demonstration of CERFAC workflows.

## File Inventory

### Primary Input Files for Workflows

#### Clinical/Case-Control Data

| File | Size | Rows | Format | Purpose | Workflow Use |
|------|------|------|--------|---------|--------------|
| **Case_control_data.tsv** | 4.0 MB | 11,228 | TSV | Case-control observations from multiple studies (BRIDGES, CARRIERS, UKB) with variant IDs in VCF format and OR calculations | INPUT: `merge_clinical_functional_data` (CLINICAL_DATA) |
| **BRCA_BRIDGES_case_control_clinical.csv** | 4.0 MB | 11,227 | CSV | BRCA2-specific case-control data from BRIDGES study, similar structure to Case_control_data | REFERENCE: BRCA-specific variant data |
| **mock_case_data.csv** | 27 KB | 1,000 | CSV | Simplified mock case data for quick testing with identifier/cdot/pdot format | TESTING: Lightweight test case |

#### Functional Assay Scores

| File | Size | Rows | Format | Purpose | Workflow Use |
|------|------|------|--------|---------|--------------|
| **Uncalibrated_functional_assay.tsv** | 366 KB | 11,010 | TSV | Functional assay scores with VCF variant IDs (chr:pos:ref:alt), amino acid changes, and raw scores | INPUT: `merge_clinical_functional_data` (FUNCTIONAL_SCORES) |
| **BRCA1_RING_BRCT_FA_sge_scores.csv** | 538 KB | 11,009 | CSV | BRCA1-specific functional assay scores from SGE (Saturation Genome Editing) experiment for RING-BRCT domain | REFERENCE: BRCA1 assay scores |
| **brca1_tavtigian_mavedb.csv** | 65 KB | 853 | CSV | BRCA1 MAVE-DB scores with HGVS notation (c. and p. format), BASE_SCALED scores, observation counts | REFERENCE: BRCA1 calibrated scores |

#### MAVEdb Export Files (Raw Data)

| File | Size | Rows | Format | Purpose |
|------|------|------|--------|---------|
| **urn_mavedb_00000003-a-1_counts.csv** | 7.9 MB | 20,725 | CSV | Raw experimental counts from MAVE-DB export with multiple replicates and conditions |
| **urn_mavedb_00000003-a-1_scores.csv** | 5.9 MB | 20,724 | CSV | Calculated scores derived from the counts file with normalized/scaled values |
| **urn_mavedb_00000006-a-1_metadata.txt** | 64 B | 1 | JSON | Metadata for MAVEDB export (chromosome 9, coordinates in hg19) |

## Data Format Details

### Clinical Data (Case_control_data.tsv)

**Key columns:**
- First column: Variant ID in VCF format `CHR-POS-REF-ALT` (e.g., `13-32314431-T-C`)
- Genomic info: CHR, POS__hg38, REF, ALT, DISTANCE, Consequence
- Study data: Multiple columns per study (BRIDGES, CARRIERS, UKB) with:
  - `_N`: Total participants
  - `_Case_Carriers_N__pct`: Carrier frequency in cases
  - `_Control_Carriers_N__pct`: Carrier frequency in controls
  - `_OR__95pct_CI`: Odds ratio with confidence interval
  - `_p-value`: Statistical significance

**Format requirements:**
- Tab-separated
- First column is variant ID (VCF format)
- Header row required

### Functional Assay Scores (Uncalibrated_functional_assay.tsv)

**Key columns:**
- First column: Variant ID in VCF format
- `chr`, `pos`, `ref`, `alt`: Genomic coordinates (also in ID)
- `aaref`, `aaalt`, `aapos`: Amino acid reference, alternative, position
- `score`: Functional assay score (may be NA)

**Format requirements:**
- Tab-separated
- First column is variant ID
- Can have NA values for missing scores

### gnomAD Variants (Expected from get_gnomad_variants workflow)

**Expected format from `get_gnomad_variants.wdl`:**
- CSV with columns: `CERFAC_variant_id_VCF`, `HGVS_cDNA_ID`, population frequencies, variant effects
- First column: VCF format variant ID
- Used as INPUT to `merge_clinical_functional_data`

## Workflow Data Flow

```
┌──────────────────────────────────┐
│  get_gnomad_variants workflow    │
│  (generates calibration variants)│
└──────────────┬───────────────────┘
               │ Output: gnomad_variants.csv
               ↓
      ┌────────────────────────────┐
      │ merge_clinical_functional  │
      │ _data workflow             │
      │ ──────────────────────────  │
      │ Inputs:                    │
      │ • gnomad_variants.csv      │
      │ • Case_control_data.tsv◄───┼── example_data/
      │ • functional_scores.tsv◄───┼── example_data/
      └────────────────────────────┘
               │ Output: merged_variants.csv
               ↓
      ┌────────────────────────────┐
      │ R Jupyter Notebook         │
      │ (Statistical calibration)  │
      └────────────────────────────┘
```

## Using Example Data for Local Testing

### Quick Test (BRCA1)

```bash
# Copy example assay scores and clinical data
cp example_data/Uncalibrated_functional_assay.tsv ./assay.tsv
cp example_data/Case_control_data.tsv ./clinical.tsv

# Run gnomAD workflow (if cloud credentials available)
java -jar tools/cromwell-92.jar run \
  workflows/get_gnomad_variants/get_gnomad_variants.wdl \
  --inputs workflows/test/get_gnomad_variants.brca1.input.json

# Run merge workflow with example data
java -jar tools/cromwell-92.jar run \
  workflows/combined_gnomad_clinvar/merge_clinical_functional_data.wdl \
  --inputs merge_test.json
```

### Test Input JSON (merge_test.json)

```json
{
  "merge_clinical_data.GENE_NAME": "BRCA1",
  "merge_clinical_data.VARIANTS_FILE": "path/to/BRCA1_gnomad_variants.csv",
  "merge_clinical_data.FUNCTIONAL_SCORES": "example_data/Uncalibrated_functional_assay.tsv",
  "merge_clinical_data.CLINICAL_DATA": "example_data/Case_control_data.tsv"
}
```

## Data Characteristics

### Variant ID Formats in Example Data

The example data uses **VCF format** variant IDs:
- Format: `CHR-POS-REF-ALT` (dashes as separators)
- Example: `17-43099873-G-A` (chromosome 17, position 43099873, G→A)
- Supported by `merge_clinical_functional_data` workflow

### Coverage

- **BRCA1/BRCA2 variants**: Examples cover both genes with real clinical data
- **gnomAD frequencies**: Included in Case_control_data.tsv
- **Case/control counts**: Real data from BRIDGES, CARRIERS, and UKB studies
- **Functional assay types**: SGE, Y2H, and MAVEdb formats included

## Licensing & Attribution

Example data sources:
- **BRIDGES study**: Case-control data for BRCA variants
- **MAVEdb**: Multiplexed assay of variant effect database
- **gnomAD**: Genome Aggregation Database (population frequencies)

For publication, cite the respective data sources.

## Notes for Workflow Developers

1. **Data must be reformatted** to match workflow input expectations
2. **Variant ID in first column** is critical for merge workflow
3. **Header row required** for all CSV/TSV files
4. **NA/missing values** are handled by merge workflow (skipped)
5. **Case/control proportions** are provided as percentages or counts depending on data source
