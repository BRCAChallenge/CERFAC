# Docker Images

Three Docker images power the CERFAC pipeline, hosted on Docker Hub at [`brcachallenge/cerfac`](https://hub.docker.com/r/brcachallenge/cerfac):

## 1. brcachallenge/cerfac:clinvar-latest
**Purpose**: Extract ClinVar variant data from NCBI  
**Base**: ubuntu:22.04  
**Key tools**: NCBI EDirect, Python (pandas, natsort)  
**Size**: ~257 MB compressed

Runs `get_clinvar_variants.wdl` to query ClinVar via EDirect and parse XML responses.

## 2. brcachallenge/cerfac:merge-latest
**Purpose**: Merge functional assay scores with clinical data  
**Base**: ubuntu:22.04  
**Key tools**: Python (pandas, natsort, requests)  
**Size**: ~242 MB compressed

Runs `merge_clinical_functional_data.wdl` to normalize variants and merge datasets via VRS API.

## 3. brcachallenge/cerfac:gnomad-latest
**Purpose**: Query population variant frequencies  
**Base**: hailgenetics/hail:0.2.131-py3.11  
**Key tools**: Hail, gnomad package, Python  
**Size**: ~1.28 GB compressed

Runs `get_gnomad_vars.wdl` to retrieve allele frequencies from gnomAD.

## Using Published Images

Images are automatically pulled from Docker Hub when running workflows on Terra:

```bash
docker pull brcachallenge/cerfac:clinvar-latest
docker pull brcachallenge/cerfac:merge-latest
docker pull brcachallenge/cerfac:gnomad-latest
```

## Building Locally

To build images locally for development:

```bash
./docker/build.sh                    # Build for current platform
./docker/build.sh linux/amd64        # Build for AMD64 (production)
./docker/build.sh linux/arm64        # Build for ARM64 (development)
```

Or manually:

```bash
docker build -t cerfac-clinvar docker/clinvar/
docker build -t cerfac-merge docker/merge_clinical_data/
docker build -t cerfac-gnomad workflows/get_gnomad_variants/docker/
```

## Testing

Validate image functionality:

```bash
docker run --rm brcachallenge/cerfac:clinvar-latest esearch -help
docker run --rm brcachallenge/cerfac:merge-latest python3 -c "import requests; print('OK')"
docker run --rm brcachallenge/cerfac:gnomad-latest python3 -c "import hail; print('OK')"
```

Or with local builds:

```bash
docker run --rm cerfac-clinvar esearch -help
docker run --rm cerfac-merge python3 -c "import requests; print('OK')"
docker run --rm cerfac-gnomad python3 -c "import hail; print('OK')"
```

## Notes

- ClinVar and Merge images are 95% identical; consolidation possible in future
- gnomAD image is large due to Hail framework
- EDirect downloads from NCBI at build time (~1-2 minutes)
- All images assume docker user created (though workflows run as root)
