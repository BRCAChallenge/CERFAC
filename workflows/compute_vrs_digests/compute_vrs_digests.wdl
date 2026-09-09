version 1.0

workflow compute_vrs_digests {
    meta {
        author: "CERFAC Team"
        description: "Compute VRS digests and merge multiple variant datasets"
    }

    input {
        String gene_name
        File gnomad_variants_file
        File? functional_scores_file
        File? clinical_data_file
        String assembly = "GRCh38"
        String docker_image = "brcachallenge/cerfac:vrs-digests-latest"
    }

    call compute_gnomad_vrs_digests {
        input:
            gene_name = gene_name,
            input_file = gnomad_variants_file,
            assembly = assembly,
            docker_image = docker_image
    }

    call compute_functional_assay_vrs_digests {
        input:
            gene_name = gene_name,
            input_file = functional_scores_file,
            assembly = assembly,
            docker_image = docker_image
    }

    call compute_clinical_vrs_digests {
        input:
            gene_name = gene_name,
            input_file = clinical_data_file,
            assembly = assembly,
            docker_image = docker_image
    }

    call merge_vrs_data {
        input:
            gene_name = gene_name,
            gnomad_file = compute_gnomad_vrs_digests.output_file,
            functional_file = compute_functional_assay_vrs_digests.output_file,
            clinical_file = compute_clinical_vrs_digests.output_file,
            docker_image = docker_image
    }

    output {
        File gnomad_with_vrs = compute_gnomad_vrs_digests.output_file
        File? functional_assay_with_vrs = compute_functional_assay_vrs_digests.output_file
        File? clinical_with_vrs = compute_clinical_vrs_digests.output_file
        File merged_data = merge_vrs_data.merged_output
    }
}


task compute_gnomad_vrs_digests {
    input {
        String gene_name
        File input_file
        String assembly = "GRCh38"
        String docker_image
        Int memory_gb = 4
        Int disk_gb = 20
    }

    String output_filename = gene_name + "_gnomad_with_vrs.csv"
    String seqrepo_path = "seqrepo+file:///seqrepo-" + assembly + "/master"

    command <<<
        set -eux -o pipefail
        python3 /vrs-python/compute_vrs_digests.py \
            --input ~{input_file} \
            --output ~{output_filename} \
            --variant-type gnomad \
            --seqrepo-uri ~{seqrepo_path} \
            --assembly ~{assembly}
    >>>

    output {
        File output_file = output_filename
    }

    runtime {
        docker: docker_image
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_gb + " SSD"
        cpu: 1
    }
}


task compute_functional_assay_vrs_digests {
    input {
        String gene_name
        File? input_file
        String assembly = "GRCh38"
        String docker_image
        Int memory_gb = 4
        Int disk_gb = 20
    }

    String output_filename = gene_name + "_functional_assay_with_vrs.csv"
    String seqrepo_path = "seqrepo+file:///seqrepo-" + assembly + "/master"

    command <<<
        set -eux -o pipefail
        if [ -z "~{input_file}" ]; then
            touch ~{output_filename}
            exit 0
        fi
        python3 /vrs-python/compute_vrs_digests.py \
            --input ~{input_file} \
            --output ~{output_filename} \
            --variant-type functional_assay \
            --seqrepo-uri ~{seqrepo_path} \
            --assembly ~{assembly}
    >>>

    output {
        File? output_file = if defined(input_file) then output_filename else None
    }

    runtime {
        docker: docker_image
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_gb + " SSD"
        cpu: 1
    }
}


task compute_clinical_vrs_digests {
    input {
        String gene_name
        File? input_file
        String assembly = "GRCh38"
        String docker_image
        Int memory_gb = 4
        Int disk_gb = 20
    }

    String output_filename = gene_name + "_clinical_with_vrs.csv"
    String seqrepo_path = "seqrepo+file:///seqrepo-" + assembly + "/master"

    command <<<
        set -eux -o pipefail
        if [ -z "~{input_file}" ]; then
            touch ~{output_filename}
            exit 0
        fi
        python3 /vrs-python/compute_vrs_digests.py \
            --input ~{input_file} \
            --output ~{output_filename} \
            --variant-type clinical \
            --seqrepo-uri ~{seqrepo_path} \
            --assembly ~{assembly}
    >>>

    output {
        File? output_file = if defined(input_file) then output_filename else None
    }

    runtime {
        docker: docker_image
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_gb + " SSD"
        cpu: 1
    }
}


task merge_vrs_data {
    input {
        String gene_name
        File gnomad_file
        File? functional_file
        File? clinical_file
        String docker_image
        Int memory_gb = 4
        Int disk_gb = 20
    }

    String output_filename = gene_name + "_merged_vrs.csv"

    command <<<
        python3 << 'PYTHON_CODE'
import pandas as pd
import sys

gene_name = "~{gene_name}"
output_file = "~{output_filename}"

# Load gnomAD data (always present)
print(f"Loading gnomAD data...")
gnomad = pd.read_csv("~{gnomad_file}")
print(f"  gnomAD: {len(gnomad)} rows, {len(gnomad.columns)} columns")

# Rename vrs_digest to have source prefix for clarity
gnomad = gnomad.rename(columns={
    'vrs_digest': 'gnomad_vrs_digest',
    'vrs_id': 'gnomad_vrs_id',
    'genomic_hgvs': 'gnomad_genomic_hgvs'
})

# Prepare for merge - use vrs_id as merge key (using the renamed gnomad version)
gnomad = gnomad.rename(columns={'gnomad_vrs_id': 'vrs_id'})

merged = gnomad.copy()

# Load and merge functional scores if provided
if "~{functional_file}" != "":
    print(f"Loading functional assay data...")
    functional = pd.read_csv("~{functional_file}")
    print(f"  Functional: {len(functional)} rows, {len(functional.columns)} columns")

    # Keep only variant_id and VRS columns for merge
    functional_cols = [col for col in functional.columns if col in ['vrs_id', 'vrs_digest', 'genomic_hgvs']]
    functional_cols_renamed = {
        'vrs_digest': 'functional_vrs_digest',
        'genomic_hgvs': 'functional_genomic_hgvs'
    }
    functional = functional[[c for c in functional.columns if 'vrs' in c.lower() or 'hgvs' in c.lower() or c == 'variant_id_vcf']].copy()
    functional = functional.rename(columns=functional_cols_renamed)

    # Merge on vrs_id
    merged = pd.merge(merged, functional, on='vrs_id', how='outer', suffixes=('_gnomad', '_functional'))
    print(f"  After merge: {len(merged)} rows")

# Load and merge clinical data if provided
if "~{clinical_file}" != "":
    print(f"Loading clinical data...")
    clinical = pd.read_csv("~{clinical_file}")
    print(f"  Clinical: {len(clinical)} rows, {len(clinical.columns)} columns")

    clinical_cols_renamed = {
        'vrs_digest': 'clinical_vrs_digest',
        'genomic_hgvs': 'clinical_genomic_hgvs'
    }
    clinical = clinical[[c for c in clinical.columns if 'vrs' in c.lower() or 'hgvs' in c.lower() or c == 'vrs_id']].copy()
    clinical = clinical.rename(columns=clinical_cols_renamed)

    # Merge on vrs_id
    merged = pd.merge(merged, clinical, on='vrs_id', how='outer', suffixes=('_gnomad', '_clinical'))
    print(f"  After merge: {len(merged)} rows")

# Save merged data
merged.to_csv(output_file, index=False)
print(f"\nMerged data saved to {output_file}")
print(f"  Total rows: {len(merged)}")
print(f"  Total columns: {len(merged.columns)}")

PYTHON_CODE
    >>>

    output {
        File merged_output = output_filename
    }

    runtime {
        docker: docker_image
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_gb + " SSD"
        cpu: 1
    }
}
