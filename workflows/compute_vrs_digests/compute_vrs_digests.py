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
import itertools
import os
import re
from pathlib import Path

# Disable bioutils network retries for offline operation
os.environ['BIOUTILS_NCBI_RETRIES'] = '0'
os.environ['HGVS_SEQREPO_DIR'] = '/seqrepo-GRCh38'

from ga4gh.vrs.dataproxy import create_dataproxy
from ga4gh.vrs.extras.translator import AlleleTranslator

VCF_ID_PATTERN = re.compile(r'^(chr)?[0-9XYM]{1,2}-[0-9]+-[ACGTN]+-[ACGTN]+$', re.IGNORECASE)
LOCUS_ALIASES = {
    'chr': ('chr', 'chrom', 'chromosome'),
    'pos': ('pos', 'position', 'start'),
    'ref': ('ref',),
    'alt': ('alt',),
}


def sniff_delimiter(path: str) -> str:
    with open(path, 'r', encoding='utf-8-sig') as f:
        return '\t' if '\t' in f.readline() else ','


def detect_variant_id_column(path: str, delimiter: str):
    """Return the column whose sampled values are mostly VCF-style IDs, or None."""
    with open(path, 'r', encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        reader.fieldnames = [name.replace('﻿', '') for name in reader.fieldnames]
        sample = list(itertools.islice(reader, 50))

    best, best_fraction = None, 0.0
    for name in reader.fieldnames:
        values = [str(row.get(name) or '').strip() for row in sample]
        values = [v for v in values if v]
        if not values:
            continue
        fraction = sum(bool(VCF_ID_PATTERN.match(v)) for v in values) / len(values)
        if fraction > best_fraction:
            best, best_fraction = name, fraction
    return best if best_fraction >= 0.9 else None


def resolve_locus_columns(fieldnames, requested=None):
    """Map chr/pos/ref/alt to column names, from an explicit list or by alias matching."""
    if requested:
        names = [name.strip() for name in requested.split(',')]
        if len(names) != 4 or not set(names).issubset(fieldnames):
            raise ValueError(f"--locus-columns must be 4 names chr,pos,ref,alt present in the file, got {requested!r}")
        return dict(zip(['chr', 'pos', 'ref', 'alt'], names))

    lowered = {name.lower(): name for name in fieldnames}
    locus = {}
    for key, aliases in LOCUS_ALIASES.items():
        match = next((lowered[alias] for alias in aliases if alias in lowered), None)
        if match is None:
            return None
        locus[key] = match
    return locus


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
        try:
            allele = translator.translate_from(variant_str, variant_format)
        except Exception as e:
            return {
                "variant_id": variant_str,
                "vrs_digest": None,
                "vrs_id": None,
                "genomic_hgvs": None,
                "vrs_start": None,
                "vrs_end": None,
                "error": f"translate_from: {type(e).__name__}"
            }

        # Extract genomic HGVS representation
        try:
            genomic_hgvs_list = translator.translate_to(allele, "hgvs")
            genomic_hgvs = genomic_hgvs_list[0] if genomic_hgvs_list else None
        except Exception as e:
            # Still return the allele even if HGVS translation fails
            genomic_hgvs = None

        # Extract location information
        try:
            location = allele.location
            pos_start = location.start if location else None
            pos_end = location.end if location else None
        except Exception as e:
            pos_start = None
            pos_end = None

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
            "error": f"outer: {type(e).__name__}: {str(e)[:50]}"
        }


def process_gnomad_variants(input_file: str, translator: AlleleTranslator, output_file: str, variant_id_column: str = None):
    """
    Process gnomAD variants with VCF format IDs.

    Input: CSV with CERFAC_variant_id_VCF_gnomad column (or variant_id_column if given)
    Output: CSV with original columns + vrs_digest, genomic_hgvs columns
    """
    print(f"Processing gnomAD variants from {input_file}")

    rows_processed = 0
    rows_with_digest = 0
    rows_with_error = 0

    with open(input_file, 'r', encoding='utf-8-sig') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        reader.fieldnames = [name.replace('﻿', '') for name in reader.fieldnames]

        id_col = variant_id_column or (
            'CERFAC_variant_id_VCF_gnomad' if 'CERFAC_variant_id_VCF_gnomad' in reader.fieldnames
            else detect_variant_id_column(input_file, ','))
        if id_col not in reader.fieldnames:
            raise ValueError(f"No variant ID column found (set --variant-id-column); got {id_col!r}")
        print(f"  Using column '{id_col}' for variant IDs")

        fieldnames = reader.fieldnames + ['genomic_hgvs', 'vrs_digest', 'vrs_id', 'vrs_start', 'vrs_end', 'vrs_error']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            variant_id = row[id_col]

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


def process_functional_assay_variants(input_file: str, translator: AlleleTranslator, output_file: str,
                                      variant_id_column: str = None, locus_columns: str = None):
    """
    Process functional assay variants.

    Variant IDs come from variant_id_column, else a detected VCF-ID column, else chr/pos/ref/alt
    (or locus_columns) columns.
    Output: CSV with original columns + variant_id_vcf, genomic_hgvs, vrs_digest columns
    """
    print(f"Processing functional assay variants from {input_file}")

    # Detect delimiter
    with open(input_file, 'r', encoding='utf-8-sig') as f:
        first_line = f.readline()
        delimiter = '\t' if '\t' in first_line else ','

    rows_processed = 0
    rows_with_digest = 0
    rows_with_error = 0

    with open(input_file, 'r', encoding='utf-8-sig') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile, delimiter=delimiter)
        reader.fieldnames = [name.replace('﻿', '') for name in reader.fieldnames]

        if variant_id_column:
            id_col = variant_id_column
        elif locus_columns is None:
            id_col = detect_variant_id_column(input_file, delimiter)
        else:
            id_col = None

        locus = None
        if id_col:
            if id_col not in reader.fieldnames:
                raise ValueError(f"Input file missing variant ID column '{id_col}'")
            print(f"  Using column '{id_col}' for variant IDs")
        else:
            locus = resolve_locus_columns(reader.fieldnames, locus_columns)
            if locus is None:
                raise ValueError("No variant ID column detected and no chr/pos/ref/alt columns found")
            print(f"  Building variant IDs from columns {locus}")

        fieldnames = reader.fieldnames + ['variant_id_vcf', 'genomic_hgvs', 'vrs_digest', 'vrs_id', 'vrs_start', 'vrs_end', 'vrs_error']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter=',')
        writer.writeheader()

        for row in reader:
            if id_col:
                variant_id_vcf = str(row[id_col]).strip()
                missing = variant_id_vcf.lower() in ('', 'na', 'none')
                missing_msg = 'Empty or NA variant'
            else:
                parts = [str(row[locus[key]]).strip() for key in ('chr', 'pos', 'ref', 'alt')]
                variant_id_vcf = f"{parts[0]}-{parts[1]}-{parts[2].upper()}-{parts[3].upper()}"
                missing = any(part.lower() in ('na', 'none', '') for part in parts)
                missing_msg = 'Missing allele information'
            row['variant_id_vcf'] = variant_id_vcf

            if missing:
                row['genomic_hgvs'] = None
                row['vrs_digest'] = None
                row['vrs_id'] = None
                row['vrs_start'] = None
                row['vrs_end'] = None
                row['vrs_error'] = missing_msg
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


def process_clinical_variants(input_file: str, translator: AlleleTranslator, output_file: str, variant_id_column: str = None):
    """
    Process clinical data variants.

    Input: TSV/CSV with variant IDs in the second column (VCF format), or variant_id_column if given
    Output: CSV with original columns + genomic_hgvs, vrs_digest columns
    """
    print(f"Processing clinical variants from {input_file}")

    # Detect delimiter
    with open(input_file, 'r', encoding='utf-8-sig') as f:
        first_line = f.readline()
        delimiter = '\t' if '\t' in first_line else ','

    rows_processed = 0
    rows_with_digest = 0
    rows_with_error = 0

    with open(input_file, 'r', encoding='utf-8-sig') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile, delimiter=delimiter)
        reader.fieldnames = [name.replace('﻿', '') for name in reader.fieldnames]

        if len(reader.fieldnames) < 2:
            raise ValueError("Input file must have at least 2 columns")

        variant_id_col = variant_id_column or detect_variant_id_column(input_file, delimiter) or reader.fieldnames[1]
        if variant_id_col not in reader.fieldnames:
            raise ValueError(f"Input file missing variant ID column '{variant_id_col}'")
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
    parser.add_argument('--variant-id-column', default=None,
                       help='Column holding the variant ID (default: auto-detected)')
    parser.add_argument('--locus-columns', default=None,
                       help='Comma-separated chr,pos,ref,alt columns for functional_assay input (default: auto-detected)')

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
            process_gnomad_variants(args.input, translator, args.output, args.variant_id_column)
        elif args.variant_type == 'functional_assay':
            process_functional_assay_variants(args.input, translator, args.output, args.variant_id_column, args.locus_columns)
        elif args.variant_type == 'clinical':
            process_clinical_variants(args.input, translator, args.output, args.variant_id_column)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
