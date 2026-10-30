#!/usr/bin/env python3
"""
Optimized VRS digest computation using ga4gh-vrs directly.
- No external network access
- Fast batch processing
- Pre-loaded translator
- Direct seqrepo access
"""

import sys
import os
import csv
import argparse
import logging
from pathlib import Path

# Set UTA_DB_URL BEFORE importing any hgvs modules
# Format: postgresql://[user:password@]host/database/schema
# The schema name MUST be included after the database name with a slash
os.environ["UTA_DB_URL"] = "postgresql://postgres:postgres@localhost/uta/uta_20241220"

from ga4gh.vrs.extras.translator import AlleleTranslator
from ga4gh.vrs.dataproxy import create_dataproxy

# Disable all network access and logging from dependencies
logging.getLogger("requests").setLevel(logging.CRITICAL)
logging.getLogger("urllib3").setLevel(logging.CRITICAL)
logging.getLogger("bioutils").setLevel(logging.CRITICAL)
os.environ["BIOUTILS_NCBI_RETRIES"] = "0"


def setup_translator(seqrepo_uri: str, assembly: str = "GRCh38") -> AlleleTranslator:
    """Initialize AlleleTranslator with local SeqRepo and UTA."""
    print(f"Initializing AlleleTranslator...", flush=True)
    print(f"  SeqRepo URI: {seqrepo_uri}", flush=True)
    print(f"  Assembly: {assembly}", flush=True)

    # Verify PostgreSQL UTA is available
    # (UTA_DB_URL was already set at module level to use local PostgreSQL)
    try:
        import psycopg2
        try:
            conn = psycopg2.connect(
                host="localhost",
                user="postgres",
                password="postgres",
                database="uta",
                connect_timeout=2
            )
            conn.close()
            print(f"  UTA: Using local PostgreSQL database (UTA_DB_URL={os.environ.get('UTA_DB_URL')})", flush=True)
        except psycopg2.OperationalError as e:
            print(f"  UTA: PostgreSQL not available, network will be used", flush=True)
    except ImportError:
        print(f"  UTA: psycopg2 not available", flush=True)

    # Create data proxy for sequence access
    data_proxy = create_dataproxy(uri=seqrepo_uri)

    # Create translator with local data
    # hgvs will use UTA_DB_URL environment variable for UTA access
    translator = AlleleTranslator(
        data_proxy=data_proxy,
        default_assembly_name=assembly
    )

    print("AlleleTranslator initialized successfully", flush=True)
    return translator


def compute_vrs_digest(translator: AlleleTranslator, variant_str: str) -> dict:
    """Compute VRS digest for a single variant."""
    try:
        # Translate from gnomad VCF format to VRS allele
        allele = translator.translate_from(variant_str, "gnomad")

        # Get genomic HGVS representation
        genomic_hgvs_list = translator.translate_to(allele, "hgvs")
        genomic_hgvs = genomic_hgvs_list[0] if genomic_hgvs_list else None

        # Extract location information
        location = allele.location
        vrs_start = location.start if location else None
        vrs_end = location.end if location else None

        return {
            "variant_id": variant_str,
            "vrs_digest": allele.digest,
            "vrs_id": allele.id,
            "genomic_hgvs": genomic_hgvs,
            "vrs_start": vrs_start,
            "vrs_end": vrs_end,
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


def process_variants(input_file: str, translator: AlleleTranslator, output_file: str, variant_type: str):
    """Process variants and write VRS digests."""
    print(f"Processing variants from {input_file}", flush=True)

    # Detect delimiter
    with open(input_file, 'r') as f:
        first_line = f.readline()
        delimiter = '\t' if '\t' in first_line else ','

    rows_processed = 0
    rows_with_digest = 0
    rows_with_error = 0

    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile, delimiter=delimiter)

        # Add VRS columns to output
        fieldnames = reader.fieldnames + ['vrs_id', 'vrs_digest', 'genomic_hgvs', 'vrs_start', 'vrs_end', 'vrs_error']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter=',')
        writer.writeheader()

        for row in reader:
            # Extract variant ID based on type
            if variant_type == 'gnomad':
                # For gnomAD: look for CERFAC_variant_id_VCF_gnomad column
                variant_id = str(row.get('CERFAC_variant_id_VCF_gnomad', row[reader.fieldnames[-1]])).strip()
            elif variant_type == 'functional_assay':
                # For functional assay: construct from chr/pos/ref/alt
                chr_val = str(row.get('chr', '')).strip()
                pos_val = str(row.get('pos', '')).strip()
                ref_val = str(row.get('ref', '')).strip().upper()
                alt_val = str(row.get('alt', '')).strip().upper()
                variant_id = f"{chr_val}-{pos_val}-{ref_val}-{alt_val}"
            elif variant_type == 'clinical':
                # For clinical: second column has variant IDs
                variant_id = str(row[reader.fieldnames[1]]).strip()
            else:
                variant_id = str(row[reader.fieldnames[-1]]).strip()

            # Compute VRS digest
            if variant_id and not variant_id.lower() in ('na', 'none', ''):
                result = compute_vrs_digest(translator, variant_id)
                row['vrs_id'] = result['vrs_id']
                row['vrs_digest'] = result['vrs_digest']
                row['genomic_hgvs'] = result['genomic_hgvs']
                row['vrs_start'] = result['vrs_start']
                row['vrs_end'] = result['vrs_end']
                row['vrs_error'] = result['error']

                if result['vrs_digest']:
                    rows_with_digest += 1
                else:
                    rows_with_error += 1
            else:
                row['vrs_id'] = None
                row['vrs_digest'] = None
                row['genomic_hgvs'] = None
                row['vrs_start'] = None
                row['vrs_end'] = None
                row['vrs_error'] = 'Empty variant ID'
                rows_with_error += 1

            writer.writerow(row)
            rows_processed += 1

            if rows_processed % 100 == 0:
                print(f"  Processed {rows_processed} variants ({rows_with_digest} with digest, {rows_with_error} errors)", flush=True)

    print(f"Complete: {rows_processed} total, {rows_with_digest} with digest, {rows_with_error} errors", flush=True)
    print(f"Output written to {output_file}", flush=True)


def main():
    parser = argparse.ArgumentParser(description="Compute VRS digests (optimized, local-only)")
    parser.add_argument('--input', required=True, help='Input file path')
    parser.add_argument('--output', required=True, help='Output file path')
    parser.add_argument('--variant-type', default='gnomad',
                       choices=['gnomad', 'functional_assay', 'clinical'],
                       help='Type of variant data')
    parser.add_argument('--seqrepo-uri', default='seqrepo+file:///seqrepo-GRCh38/master',
                       help='SeqRepo URI (must be local file:// only)')
    parser.add_argument('--assembly', default='GRCh38',
                       help='Reference assembly')

    args = parser.parse_args()

    # Verify seqrepo is local
    if not args.seqrepo_uri.startswith('seqrepo+file://'):
        raise ValueError(f"SeqRepo URI must use file:// scheme, got: {args.seqrepo_uri}")

    # Setup translator
    translator = setup_translator(args.seqrepo_uri, args.assembly)

    # Process variants
    process_variants(args.input, translator, args.output, args.variant_type)


if __name__ == '__main__':
    main()
