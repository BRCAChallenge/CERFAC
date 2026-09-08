# Cromwell Setup and Testing Guide

## Installation

Pre-downloaded Cromwell JARs are included in this repository:
- **Location**: `tools/cromwell-92.jar` (recommended) and `tools/cromwell-86.jar` (legacy)
- **Size**: ~250+ MB
- **Current Version**: Cromwell 92 (tested and verified with offline VRS workflows)

### Prerequisites

1. **Java Runtime Environment (JRE)** - Required to run Cromwell
   ```bash
   # Check if Java is installed
   java -version
   
   # macOS (using Homebrew)
   brew install java
   
   # Linux (Ubuntu/Debian)
   sudo apt-get install default-jre
   
   # Or install OpenJDK
   brew install openjdk
   ```

2. **Docker** - Required for containerized workflow execution
   ```bash
   # Check if Docker is installed
   docker --version
   ```

## Configuration

The default Cromwell configuration is in `workflows/combined_gnomad_clinvar/cromwell.conf`:

```conf
backend {
  default = Local
  providers {
    Local {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        run-in-background = true
        submit-docker = """
        docker run \
          --cidfile ${docker_cid} \
          -i \
          --memory=7g \
          --cpus=2 \
          --entrypoint ${job_shell} \
          -v ${cwd}:${docker_cwd} \
          ${docker} ${docker_script}
        """
      }
    }
  }
}
```

**Configuration Details**:
- **Backend**: Local (runs on the machine, not cloud)
- **Memory Default**: 7 GB per task
- **CPU Default**: 2 cores per task
- **Docker**: Enabled with volume mounts
- **Docker Image**: Uses `cerfac:vrs-offline` for VRS tasks (completely offline, no network dependencies)
- **Volume Mounts**: Maps current working directory and `/tmp` for data I/O

**For Offline VRS Workflows**:
The configuration automatically uses the `cerfac:vrs-offline` Docker image for VRS computation tasks, which includes:
- SeqRepo GRCh38 (local sequence reference)
- UTA database (PostgreSQL)
- ga4gh-vrs with all dependencies
- bioutils network retries disabled
- Zero external network dependencies

## Running Workflows

### Syntax Validation

Validate a WDL file without running it:

```bash
java -jar tools/cromwell-92.jar validate workflows/combined_gnomad_clinvar/merge_clinical_functional_data_vrs_offline.wdl
```

### Running a Workflow

Execute a workflow with input JSON:

```bash
cd /path/to/CERFAC
java -Dconfig.file=workflows/combined_gnomad_clinvar/cromwell.conf \
     -jar tools/cromwell-92.jar run \
     workflows/combined_gnomad_clinvar/merge_clinical_functional_data_vrs_offline.wdl \
     --inputs <input-file.json>
```

### Offline VRS Workflow (Production)

The `merge_clinical_functional_data_vrs_offline.wdl` workflow computes VRS digests **without any network access**:

```bash
# Create test inputs
cat > /tmp/test_wdl_input.json << 'EOF'
{
  "merge_clinical_data.GENE_NAME": "BRCA1",
  "merge_clinical_data.VARIANTS_FILE": "/path/to/variants.csv",
  "merge_clinical_data.FUNCTIONAL_SCORES": "/path/to/functional.csv",
  "merge_clinical_data.CLINICAL_DATA": "/path/to/clinical.csv"
}
EOF

# Run offline VRS workflow
java -Dconfig.file=workflows/combined_gnomad_clinvar/cromwell.conf \
     -jar tools/cromwell-92.jar run \
     workflows/combined_gnomad_clinvar/merge_clinical_functional_data_vrs_offline.wdl \
     --inputs /tmp/test_wdl_input.json
```

**Workflow Tasks**:
- `vrs_variants` - Process gnomAD format variants (chr-pos-ref-alt)
- `vrs_scores` - Process functional assay variants (separate chr/pos/ref/alt columns)
- `vrs_clinical` - Process clinical data variants
- `merge_vrs_files` - Merge all VRS-annotated files into single output

**Expected Output**:
```
Workflow 4a07d53c-b884-40ba-b2b1-1f65b6bd6983 transitioned to state Succeeded
```

**Performance**:
- Per-variant: 0.3-0.5 seconds
- 26,000 variants: ~2-3 hours
- Zero network calls (completely offline)

## Common Commands

```bash
# Validate WDL syntax
java -jar tools/cromwell-92.jar validate workflows/combined_gnomad_clinvar/merge_clinical_functional_data_vrs_offline.wdl

# Run offline VRS workflow
java -Dconfig.file=workflows/combined_gnomad_clinvar/cromwell.conf \
     -jar tools/cromwell-92.jar run \
     workflows/combined_gnomad_clinvar/merge_clinical_functional_data_vrs_offline.wdl \
     --inputs <input-file.json>

# List version
java -jar tools/cromwell-92.jar --version

# Get help
java -jar tools/cromwell-92.jar --help

# Check Cromwell executions (outputs and logs)
ls -la cromwell-executions/merge_clinical_data/
tail -100 cromwell-executions/merge_clinical_data/<workflow-id>/call-vrs_variants/execution/stdout
```

## Offline VRS Computation

The `cerfac:vrs-offline` Docker image enables completely offline variant analysis:

### Prerequisites

Build the offline VRS Docker image (includes 13.3 GB of data):

```bash
cd workflows/compute_vrs_digests
docker build -t cerfac:vrs-offline .
```

**Build Time**: ~10-15 minutes  
**Image Size**: 13.3 GB (includes SeqRepo GRCh38 + UTA database + PostgreSQL)

### Features

- **Zero Network Calls**: Completely offline operation, no NCBI/biocommons lookups
- **PostgreSQL UTA**: Embedded transcript database
- **SeqRepo GRCh38**: Local sequence reference (800MB+ of genomic data)
- **Multi-format Support**: gnomAD, functional assay, clinical data formats
- **Performance**: 0.3-0.5 seconds per variant (26K variants in 2-3 hours)

### Environment Variables (Auto-configured)

```bash
GA4GH_VRS_DATAPROXY_URI=seqrepo+file:///seqrepo-GRCh38/master
UTA_DB_URL=postgresql://postgres:postgres@localhost/uta/uta_20241220
BIOUTILS_NCBI_RETRIES=0
```

### Example: Full Workflow Execution

```bash
# Test data (replace with real variant files)
cat > /tmp/variants.json << 'EOF'
{
  "merge_clinical_data.GENE_NAME": "BRCA1",
  "merge_clinical_data.VARIANTS_FILE": "/tmp/gnomad_variants.csv",
  "merge_clinical_data.FUNCTIONAL_SCORES": "/tmp/functional_scores.csv",
  "merge_clinical_data.CLINICAL_DATA": "/tmp/clinical_data.csv"
}
EOF

# Run offline workflow
java -Dconfig.file=workflows/combined_gnomad_clinvar/cromwell.conf \
     -jar tools/cromwell-92.jar run \
     workflows/combined_gnomad_clinvar/merge_clinical_functional_data_vrs_offline.wdl \
     --inputs /tmp/variants.json

# Check results
ls cromwell-executions/merge_clinical_data/*/call-merge_vrs_files/execution/
cat cromwell-executions/merge_clinical_data/*/call-merge_vrs_files/execution/stdout
```

## Troubleshooting

### Java Not Found
```
The operation couldn't be completed. Unable to locate a Java Runtime.
```
**Solution**: Install JRE/JDK (see Prerequisites above)

### Docker Not Available
```
ERROR: error creating mount: mkdir /var/lib/docker/volumes/...
```
**Solution**: Ensure Docker daemon is running and user has permissions

### OOM (Out of Memory)
Edit `cromwell.conf` to increase `--memory` in the `submit-docker` section.

### Task Timeout
Add timeout configuration to the WDL task:
```wdl
runtime {
  memory: "8 GB"
  cpu: 2
  docker: "image:tag"
  maxRetries: 3
  timeout: 3600  # seconds
}
```

## Output

Cromwell creates a `cromwell-executions/` directory containing:
- Logs: `workflow.log`, task logs
- Outputs: Task and workflow outputs
- Metadata: `metadata.json` with full execution details

## Useful Links

- [Cromwell Documentation](https://cromwell.readthedocs.io/)
- [WDL Spec](https://github.com/openwdl/wdl)
- [Cromwell GitHub Releases](https://github.com/broadinstitute/cromwell/releases)

