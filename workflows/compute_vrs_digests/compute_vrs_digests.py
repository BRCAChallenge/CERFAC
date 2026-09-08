#!/usr/bin/env python3
"""
Compute VRS digests for variants from multiple input formats.

Handles:
- gnomAD VCF format (chr-pos-ref-alt)
- Functional assay format (separate chr/pos/ref/alt columns)
- Clinical data format (variant IDs in second column)
"""

import sys
import csv
import argparse
import os
from pathlib import Path

# Disable bioutils network retries for offline operation
os.environ['BIOUTILS_NCBI_RETRIES'] = '0'
os.environ['HGVS_SEQREPO_DIR'] = '/seqrepo-GRCh38'

from ga4gh.vrs.dataproxy import create_dataproxy
from ga4gh.vrs.extras.translator import AlleleTranslator


def setup_translator(seqrepo_uri: str, assembly: str = "GRCh38") -> AlleleTranslator:
    """Initialize AlleleTranslator with SeqRepo and UTA."""
    # Create data proxy for sequence access
    data_proxy = create_dataproxy(uri=seqrepo_uri)

    # Create translator with specified assembly
    translator = AlleleTranslator(
        data_proxy=data_proxy,
        default_assembly_name=assembly
    )

    return translator


def compute_vrs_digest(translator: AlleleTranslator, variant_str: str, variant_format: str) -> dict:
    """
    Compute VRS digest and genomic HGVS for a single variant.

    Args:
        translator: AlleleTranslator instance
        variant_str: Variant string in the specified format
        variant_format: Format type ('gnomad', 'hgvs', 'spdi', 'beacon')

    Returns:
        Dictionary with variant_id, vrs_digest, genomic_hgvs, and location info
    """
    try:
        allele = translator.translate_from(variant_str, variant_format)

        # Extract genomic HGVS representation
        genomic_hgvs_list = translator.translate_to(allele, "hgvs")
        genomic_hgvs = genomic_hgvs_list[0] if genomic_hgvs_list else None

        # Extract location information
        location = allele.location
        chrom = None
        pos_start = None
        pos_end = None

        if location and location.sequenceReference:
            # Try to extract chromosome from sequence reference
            # This is approximate - full resolution requires sequence mapping
            seq_ref = location.sequenceReference.refgetAccession
            pos_start = location.start
            pos_end = location.end

        return {
            "variant_id": variant_str,
            "vrs_digest": allele.digest,
            "vrs_id": allele.id,
            "genomic_hgvs": genomic_hgvs,
            "vrs_start": pos_start,
            "vrs_end": pos_end,
            "error": None
        }
    except Exception as e:
        return {
            "variant_id": variant_str,
            "vrs_digest": None,
            "vrs_id": None,
            "genomic_hgvs": None,
            "vrs_start": None,
            "vrs_end": None,
            "error": str(e)
        }


def process_gnomad_variants(input_file: str, translator: AlleleTranslator, output_file: str):
    """
    Process gnomAD variants with VCF format IDs.

    Input: CSV with CERFAC_variant_id_VCF_gnomad column
    Output: CSV with original columns + vrs_digest, genomic_hgvs columns
    """
    print(f"Processing gnomAD variants from {input_file}")

    rows_processed = 0
    rows_with_digest = 0
    rows_with_error = 0

    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)

        # Ensure output file has VCF column
        if 'CERFAC_variant_id_VCF_gnomad' not in reader.fieldnames:
            raise ValueError("Input file missing CERFAC_variant_id_VCF_gnomad column")

        fieldnames = reader.fieldnames + ['genomic_hgvs', 'vrs_digest', 'vrs_id', 'vrs_start', 'vrs_end', 'vrs_error']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            variant_id = row['CERFAC_variant_id_VCF_gnomad']

            # Skip empty or invalid variants
            if not variant_id or variant_id.startswith('NA') or variant_id.lower() == 'none':
                row['genomic_hgvs'] = None
                row['vrs_digest'] = None
                row['vrs_id'] = None
                row['vrs_start'] = None
                row['vrs_end'] = None
                row['vrs_error'] = 'Empty or NA variant'
                rows_with_error += 1
            else:
                result = compute_vrs_digest(translator, variant_id, 'gnomad')
                row['genomic_hgvs'] = result['genomic_hgvs']
                row['vrs_digest'] = result['vrs_digest']
                row['vrs_id'] = result['vrs_id']
                row['vrs_start'] = result['vrs_start']
                row['vrs_end'] = result['vrs_end']
                row['vrs_error'] = result['error']

                if result['vrs_digest']:
                    rows_with_digest += 1
                else:
                    rows_with_error += 1

            writer.writerow(row)
            rows_processed += 1

            if rows_processed % 100 == 0:
                print(f"  Processed {rows_processed} variants ({rows_with_digest} with digest, {rows_with_error} errors)")

    print(f"Complete: {rows_processed} total, {rows_with_digest} with digest, {rows_with_error} errors")
    print(f"Output written to {output_file}")


def process_functional_assay_variants(input_file: str, translator: AlleleTranslator, output_file: str):
    """
    Process functional assay variants with separate chr/pos/ref/alt columns.

    Input: TSV/CSV with chr, pos, ref, alt columns
    Output: CSV with original columns + variant_id_vcf, genomic_hgvs, vrs_digest columns
    """
    print(f"Processing functional assay variants from {input_file}")

    # Detect delimiter
    with open(input_file, 'r') as f:
        first_line = f.readline()
        delimiter = '\t' if '\t' in first_line else ','

    rows_processed = 0
    rows_with_digest = 0
    rows_with_error = 0

    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile, delimiter=delimiter)

        # Ensure required columns exist
        required = {'chr', 'pos', 'ref', 'alt'}
        if not required.issubset(set(reader.fieldnames)):
            raise ValueError(f"Input file missing required columns: {required}")

        fieldnames = reader.fieldnames + ['variant_id_vcf', 'genomic_hgvs', 'vrs_digest', 'vrs_id', 'vrs_start', 'vrs_end', 'vrs_error']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter=',')
        writer.writeheader()

        for row in reader:
            chr_val = str(row['chr']).strip()
            pos_val = str(row['pos']).strip()
            ref_val = str(row['ref']).strip().upper()
            alt_val = str(row['alt']).strip().upper()

            # Construct VCF format variant ID
            variant_id_vcf = f"{chr_val}-{pos_val}-{ref_val}-{alt_val}"
            row['variant_id_vcf'] = variant_id_vcf

            # Skip if any component is missing or NA
            if any(v.lower() in ('na', 'none', '') for v in [chr_val, pos_val, ref_val, alt_val]):
                row['genomic_hgvs'] = None
                row['vrs_digest'] = None
                row['vrs_id'] = None
                row['vrs_start'] = None
                row['vrs_end'] = None
                row['vrs_error'] = 'Missing allele information'
                rows_with_error += 1
            else:
                result = compute_vrs_digest(translator, variant_id_vcf, 'gnomad')
                row['genomic_hgvs'] = result['genomic_hgvs']
                row['vrs_digest'] = result['vrs_digest']
                row['vrs_id'] = result['vrs_id']
                row['vrs_start'] = result['vrs_start']
                row['vrs_end'] = result['vrs_end']
                row['vrs_error'] = result['error']

                if result['vrs_digest']:
                    rows_with_digest += 1
                else:
                    rows_with_error += 1

            writer.writerow(row)
            rows_processed += 1

            if rows_processed % 100 == 0:
                print(f"  Processed {rows_processed} variants ({rows_with_digest} with digest, {rows_with_error} errors)")

    print(f"Complete: {rows_processed} total, {rows_with_digest} with digest, {rows_with_error} errors")
    print(f"Output written to {output_file}")


def process_clinical_variants(input_file: str, translator: AlleleTranslator, output_file: str):
    """
    Process clinical data variants with variant IDs in second column.

    Input: TSV/CSV with variant IDs in second column (VCF format)
    Output: CSV with original columns + genomic_hgvs, vrs_digest columns
    """
    print(f"Processing clinical variants from {input_file}")

    # Detect delimiter
    with open(input_file, 'r') as f:
        first_line = f.readline()
        delimiter = '\t' if '\t' in first_line else ','

    rows_processed = 0
    rows_with_digest = 0
    rows_with_error = 0

    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile, delimiter=delimiter)

        # Get second column name (should be variant ID column like 'hg38')
        if len(reader.fieldnames) < 2:
            raise ValueError("Input file must have at least 2 columns")

        variant_id_col = reader.fieldnames[1]
        print(f"  Using column '{variant_id_col}' for variant IDs")

        fieldnames = reader.fieldnames + ['genomic_hgvs', 'vrs_digest', 'vrs_id', 'vrs_start', 'vrs_end', 'vrs_error']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter=',')
        writer.writeheader()

        for row in reader:
            variant_id = str(row[variant_id_col]).strip()

            # Skip empty or invalid variants
            if not variant_id or variant_id.startswith('NA') or variant_id.lower() == 'none':
                row['genomic_hgvs'] = None
                row['vrs_digest'] = None
                row['vrs_id'] = None
                row['vrs_start'] = None
                row['vrs_end'] = None
                row['vrs_error'] = 'Empty or NA variant'
                rows_with_error += 1
            else:
                result = compute_vrs_digest(translator, variant_id, 'gnomad')
                row['genomic_hgvs'] = result['genomic_hgvs']
                row['vrs_digest'] = result['vrs_digest']
                row['vrs_id'] = result['vrs_id']
                row['vrs_start'] = result['vrs_start']
                row['vrs_end'] = result['vrs_end']
                row['vrs_error'] = result['error']

                if result['vrs_digest']:
                    rows_with_digest += 1
                else:
                    rows_with_error += 1

            writer.writerow(row)
            rows_processed += 1

            if rows_processed % 100 == 0:
                print(f"  Processed {rows_processed} variants ({rows_with_digest} with digest, {rows_with_error} errors)")

    print(f"Complete: {rows_processed} total, {rows_with_digest} with digest, {rows_with_error} errors")
    print(f"Output written to {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Compute VRS digests for genomic variants")
    parser.add_argument('--input', required=True, help='Input file path')
    parser.add_argument('--output', required=True, help='Output file path')
    parser.add_argument('--variant-type', required=True,
                       choices=['gnomad', 'functional_assay', 'clinical'],
                       help='Type of variant data')
    parser.add_argument('--seqrepo-uri', default=os.environ.get('GA4GH_VRS_DATAPROXY_URI', 'refget:https://www.ncbi.nlm.nih.gov/grc/human/'),
                       help='SeqRepo URI (default: from GA4GH_VRS_DATAPROXY_URI env var or RefGet public API)')
    parser.add_argument('--assembly', default='GRCh38',
                       help='Reference assembly (default: GRCh38)')

    args = parser.parse_args()

    # Verify input file exists
    if not Path(args.input).exists():
        print(f"ERROR: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    print(f"Setting up VRS translator...")
    print(f"  SeqRepo URI: {args.seqrepo_uri}")
    print(f"  Assembly: {args.assembly}")

    try:
        translator = setup_translator(args.seqrepo_uri, args.assembly)
        print("Translator initialized successfully")
    except Exception as e:
        print(f"ERROR: Failed to initialize translator: {e}", file=sys.stderr)
        sys.exit(1)

    print()

    try:
        if args.variant_type == 'gnomad':
            process_gnomad_variants(args.input, translator, args.output)
        elif args.variant_type == 'functional_assay':
            process_functional_assay_variants(args.input, translator, args.output)
        elif args.variant_type == 'clinical':
            process_clinical_variants(args.input, translator, args.output)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
