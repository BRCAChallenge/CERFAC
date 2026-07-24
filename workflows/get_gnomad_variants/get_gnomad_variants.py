#!/usr/bin/env python3
import re
import sys

import pandas as pd
import gnomad
import hail as hl
from gnomad.resources.grch38.gnomad import public_release


# Population label -> freq_index_dict key
POPULATION_FREQ_KEYS = {
    'afr': 'afr_adj',
    'amr': 'amr_adj',
    'eas': 'eas_adj',
    'nfe': 'nfe_adj',
    'rmi': 'remaining_adj',
    'sas': 'sas_adj',
    'mid': 'mid_adj',
    'fin': 'fin_adj',
    'asj': 'asj_adj',
}

TRANSCRIPT_COLS = [
    'txpt_amino_acids', 'txpt_appris', 'txpt_biotype', 'txpt_canonical',
    'variant_effect', 'txpt_distance', 'txpt_domains', 'txpt_exon',
    'txpt_hgvsc', 'txpt_hgvsp',
    'txpt_gene_pheno', 'txpt_gene_symbol', 'txpt_impact', 'txpt_intron',
    'txpt_lof', 'txpt_lof_filter', 'txpt_lof_flags', 'txpt_lof_info',
    'txpt_mane_plus_clinical', 'txpt_mane_select', 'txpt_protein_end',
    'txpt_protein_start', 'txpt_transcript_id', 'txpt_uniprot_isoform',
]

NON_CODING_EFFECTS = {
    'upstream_gene_variant',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'intron_variant',
}

DROP_COLS = [
    'txpt_appris', 'txpt_distance', 'txpt_biotype', 'txpt_canonical',
    'txpt_gene_pheno', 'txpt_gene_symbol', 'seq_region_name', 'chr_id',
    'ref_genome', 'txpt_protein_end', 'txpt_protein_start',
]


# ---------------------------------------------------------------------------
# Hail expression helpers
# ---------------------------------------------------------------------------

def get_ref_genome(locus: hl.expr.LocusExpression) -> hl.expr.StringExpression:
    return hl.str(gnomad.utils.reference_genome.get_reference_genome(locus).name)

def normalized_contig(contig: hl.expr.StringExpression) -> hl.expr.StringExpression:
    return hl.rbind(hl.str(contig).replace("^chr", ""), lambda c: hl.if_else(c == "MT", "M", c))

def chromosome_id(locus: hl.expr.LocusExpression, alleles: hl.expr.ArrayExpression,
                  max_length: int = None) -> hl.expr.StringExpression:
    chr_id = hl.str(locus.contig)
    return chr_id[0:max_length] if max_length is not None else chr_id

def get_alt_allele(locus: hl.expr.LocusExpression, alleles: hl.expr.ArrayExpression) -> hl.expr.StringExpression:
    return hl.str(alleles[1])

def get_ref_allele(locus: hl.expr.LocusExpression, alleles: hl.expr.ArrayExpression) -> hl.expr.StringExpression:
    return hl.str(alleles[0])

def get_alt_allele_len(locus: hl.expr.LocusExpression, alleles: hl.expr.ArrayExpression) -> hl.expr.Int32Expression:
    return hl.len(alleles[1])

def get_ref_allele_len(locus: hl.expr.LocusExpression, alleles: hl.expr.ArrayExpression) -> hl.expr.Int32Expression:
    return hl.len(alleles[0])

def get_end_pos_alt(locus: hl.expr.LocusExpression, alleles: hl.expr.ArrayExpression) -> hl.expr.Int32Expression:
    return hl.int32(locus.position) + hl.int32(hl.len(alleles[1]) - 1)

def get_VRS_start_ref(pos_VRS_starts: hl.expr.ArrayNumericExpression) -> hl.expr.StringExpression:
    return hl.str(pos_VRS_starts[0])

def get_VRS_start_alt(locus: hl.expr.LocusExpression,
                      pos_VRS_starts: hl.expr.ArrayNumericExpression) -> hl.expr.StringExpression:
    return hl.str(pos_VRS_starts[1])

def get_VRS_stop_ref(pos_VRS_stops: hl.expr.ArrayNumericExpression) -> hl.expr.StringExpression:
    return hl.str(pos_VRS_stops[0])

def get_VRS_stop_alt(pos_VRS_stops: hl.expr.ArrayNumericExpression) -> hl.expr.StringExpression:
    return hl.str(pos_VRS_stops[1])

def variant_idCERFAC_VCF(locus: hl.expr.LocusExpression, alleles: hl.expr.ArrayExpression,
                          max_length: int = None) -> hl.expr.StringExpression:
    var_id = hl.str(locus.contig) + "-" + hl.str(locus.position) + "-" + alleles[0] + "-" + alleles[1]
    return var_id[0:max_length] if max_length is not None else var_id


# ---------------------------------------------------------------------------
# Hail pipeline stages
# ---------------------------------------------------------------------------

def select_fields(ht: hl.Table, n_alt_alleles_key: str) -> hl.Table:
    """Select the standard variant and transcript fields from an exomes or genomes table.

    Args:
        ht: gnomAD Hail table (exomes or genomes).
        n_alt_alleles_key: Output column name for the alt allele count
            ('n_alt_alleles_exomes' or 'n_alt_alleles_genomes').
    """
    tc = ht.vep.transcript_consequences
    return ht.select(
        ht.freq,
        ht.vep.allele_string,
        ht.vep.start,
        ht.vep.end,
        ht.vep.seq_region_name,
        ht.vep.variant_class,
        ht.in_silico_predictors.spliceai_ds_max,
        **{n_alt_alleles_key: ht.allele_info.n_alt_alleles},
        VRS_Allele_IDs=ht.info.vrs.VRS_Allele_IDs,
        allele_list_VRS=ht.info.vrs.VRS_States,
        pos_VRS_starts=ht.info.vrs.VRS_Starts,
        pos_VRS_stops=ht.info.vrs.VRS_Ends,
        txpt_amino_acids=tc.amino_acids,
        txpt_appris=tc.appris,
        txpt_biotype=tc.biotype,
        txpt_canonical=tc.canonical,
        txpt_distance=tc.distance,
        txpt_domains=tc.domains,
        txpt_exon=tc.exon,
        txpt_gene_pheno=tc.gene_pheno,
        txpt_gene_symbol=tc.gene_symbol,
        txpt_hgvsc=tc.hgvsc,
        txpt_hgvsp=tc.hgvsp,
        txpt_impact=tc.impact,
        txpt_intron=tc.intron,
        txpt_lof=tc.lof,
        txpt_lof_flags=tc.lof_flags,
        txpt_lof_filter=tc.lof_filter,
        txpt_lof_info=tc.lof_info,
        txpt_mane_plus_clinical=tc.mane_plus_clinical,
        txpt_mane_select=tc.mane_select,
        txpt_protein_end=tc.protein_end,
        txpt_protein_start=tc.protein_start,
        txpt_transcript_id=tc.transcript_id,
        txpt_uniprot_isoform=tc.uniprot_isoform,
        variant_effect=tc.consequence_terms,
    )


def annotate_pop_frequencies(ht: hl.Table, col_prefix: str) -> hl.Table:
    """Add per-population adjusted allele frequency columns.

    Args:
        ht: Hail table with freq and freq_index_dict fields.
        col_prefix: Prefix for the output columns
            (e.g. 'exome_freq_main_adj_' or 'genome_freq_adj_').
    """
    return ht.annotate(**{
        f"{col_prefix}{pop}": ht.freq[ht.freq_index_dict[freq_key]]
        for pop, freq_key in POPULATION_FREQ_KEYS.items()
    })


def filter_to_gene_region(ht: hl.Table, chr_string: str,
                           start_pos: hl.expr.Int32Expression,
                           stop_pos: hl.expr.Int32Expression) -> hl.Table:
    """Filter a Hail table to variants within the specified chromosomal interval."""
    ht = ht.filter(ht.locus.contig == chr_string)
    return hl.filter_intervals(
        ht, [hl.locus_interval(chr_string, start_pos, stop_pos, reference_genome='GRCh38')]
    )


def merge_exomes_genomes(exomes: hl.Table, genomes: hl.Table) -> hl.Table:
    """Combine exome and genome tables via inner join and anti-joins.

    Variants present in both datasets are labelled 'both'; those found only
    in one source are labelled 'exomes' or 'genomes'. Returns a single
    flattened, unioned table with an 'exorgen' provenance column.
    """
    inner = genomes.join(exomes, how='inner')
    inner = inner.annotate(exorgen="both")
    inner = inner.drop(*(x for x in inner.row if re.search(r'_1', x)))

    genomes_only = genomes.anti_join(exomes).annotate(exorgen="genomes")
    exomes_only = exomes.anti_join(genomes).annotate(exorgen="exomes")

    inner = inner.drop('freq').flatten()
    exomes_only = exomes_only.drop('freq').flatten()
    genomes_only = genomes_only.drop('freq').flatten()

    union = exomes_only.union(inner, unify=True)
    union.count()
    return union.union(genomes_only, unify=True)


def annotate_variant_ids(ht: hl.Table) -> hl.Table:
    """Add CERFAC variant ID columns and allele/VRS position fields."""
    return ht.annotate(
        ref_genome=get_ref_genome(ht.locus),
        chr_id=chromosome_id(ht.locus, ht.alleles),
        CERFAC_variant_id_VCF=variant_idCERFAC_VCF(ht.locus, ht.alleles),
        allele_ref=get_ref_allele(ht.locus, ht.alleles),
        allele_alt=get_alt_allele(ht.locus, ht.alleles),
        variant_length_ref=get_ref_allele_len(ht.locus, ht.alleles),
        variant_length_alt=get_alt_allele_len(ht.locus, ht.alleles),
        pos_start_alt_vrs=get_VRS_start_alt(ht.locus, ht.pos_VRS_starts),
        pos_start_ref_vrs=get_VRS_start_ref(ht.pos_VRS_starts),
        pos_stop_alt_vrs=get_VRS_stop_alt(ht.pos_VRS_stops),
        pos_stop_ref_vrs=get_VRS_stop_ref(ht.pos_VRS_stops),
    )


# ---------------------------------------------------------------------------
# Pandas processing
# ---------------------------------------------------------------------------

def process_dataframe(df: pd.DataFrame, gene_name: str) -> pd.DataFrame:
    """Apply transcript filtering, cleaning, and column transformations.

    Args:
        df: Raw pandas DataFrame from Hail to_pandas().
        gene_name: Gene symbol to retain (e.g. 'BRCA1').

    Returns:
        Processed DataFrame ready for output.
    """
    df = df.sort_index(axis=1)

    # Explode per-transcript list columns to one row per transcript
    df = df.explode(TRANSCRIPT_COLS)

    # Filter to canonical MANE-select transcripts for the target gene
    df = df[df.txpt_canonical == 1]
    df = df.dropna(subset=['txpt_gene_symbol'])
    df = df[df.txpt_gene_symbol == gene_name]
    df = df[df['txpt_mane_select'].str.startswith('NM')]

    # Clean variant_effect (stored as a stringified list after explode)
    df['variant_effect'] = (df['variant_effect'].astype(str)
                            .str.replace('[', '', regex=False)
                            .str.replace(']', '', regex=False)
                            .str.replace("'", '', regex=False))

    # Explode uniprot isoforms to one row each
    df = df.explode(['txpt_uniprot_isoform'])

    # Drop non-coding variants
    df = df[~df.variant_effect.isin(NON_CODING_EFFECTS)]

    # Build HGVS identifier columns
    df['txpt_hgvsc'] = df['txpt_hgvsc'].astype(str)
    df['CERFAC_variant_id_HGVS_long'] = (
        df['ref_genome'].astype(str) + ':' +
        df['locus'].astype(str) + ':' +
        df['txpt_hgvsc'].astype(str)
    )
    df['txpt_hgvsc_short'] = df['txpt_hgvsc'].str.split(pat=":", n=1, regex=False).str.get(1)
    df['hgvs_pro'] = df['txpt_hgvsp'].str.split(pat=":", n=1, regex=False).str.get(1)
    df['CERFAC_variant_id_HGVS_short'] = (
        df['ref_genome'].astype(str) + ':' +
        df['locus'].astype(str) + ':' +
        df['txpt_hgvsc_short'].astype(str)
    )
    df['CERFAC_variant_id_VCF'] = df['CERFAC_variant_id_VCF'].astype(str).replace("chr", "")

    df = df.drop(columns=DROP_COLS)
    df[['txpt_exon', 'txpt_intron']] = df[['txpt_exon', 'txpt_intron']].astype(str)

    # Rename, suffix, and label source
    df = df.rename(columns={
        "alleles": "allele_list", "locus": "pos_VCF", "exorgen": "set",
        "start": "pos_start_vep", "end": "pos_stop_vep",
    }, errors='raise')
    df = df.sort_index(axis=1)
    df['variant_source'] = "gnomAD"
    df = df.add_suffix('_gnomad')
    df['HGVS_cDNA_ID_gnomad'] = (df[['txpt_mane_select_gnomad', 'txpt_hgvsc_short_gnomad']]
                                  .astype(str).agg(':'.join, axis=1))
    return df


def write_outputs(df: pd.DataFrame, gene_name: str) -> None:
    """Write variant count to gnomadcount.txt and full table to CSV."""
    import os
    count = str(df['txpt_hgvsc_short_gnomad'].nunique())

    count_file = os.path.abspath("gnomadcount.txt")
    csv_file = os.path.abspath(f"{gene_name}_gnomad_variants_MANE.csv")

    print(f"Writing outputs to: {count_file}", flush=True)
    with open(count_file, 'w') as f:
        f.write(count)

    print(f"Writing CSV to: {csv_file}", flush=True)
    df.to_csv(csv_file, sep=',', index=False)
    print(f"Outputs written successfully", flush=True)


# ---------------------------------------------------------------------------
# Top-level orchestration
# ---------------------------------------------------------------------------

def get_gnomad_variants(gene_name: str, chr_id: str,
                        start_locus: int, end_locus: int, hail_mem_gb: int) -> None:
    import os
    try:
        if 'GOOGLE_APPLICATION_CREDENTIALS' not in os.environ:
            default_creds = os.path.expanduser('~/.config/gcloud/application_default_credentials.json')
            if os.path.exists(default_creds):
                os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = default_creds

        hl.init(spark_conf={
            'spark.driver.memory': f"{hail_mem_gb}g",
            'spark.hadoop.google.cloud.auth.type': 'APPLICATION_DEFAULT'
        })

        chr_string = f"chr{chr_id}"
        start_pos = hl.int32(start_locus)
        stop_pos = hl.int32(end_locus)

        print(f"Loading gnomAD data for {gene_name}: chr{chr_id}:{start_locus}-{end_locus}", flush=True)

        # Load, select fields, and count
        v4exomes = select_fields(public_release("exomes").ht(), 'n_alt_alleles_exomes')
        v4exomes.count()
        v4genomes = select_fields(public_release("genomes").ht(), 'n_alt_alleles_genomes')
        v4genomes.count()

        # Annotate population frequencies
        v4exomes = annotate_pop_frequencies(v4exomes, 'exome_freq_main_adj_')
        v4genomes = annotate_pop_frequencies(v4genomes, 'genome_freq_adj_')

        # Filter to gene region, merge, and annotate variant IDs
        exomes_region = filter_to_gene_region(v4exomes, chr_string, start_pos, stop_pos)
        genomes_region = filter_to_gene_region(v4genomes, chr_string, start_pos, stop_pos)
        gnomad_union = merge_exomes_genomes(exomes_region, genomes_region)
        gnomad_union = annotate_variant_ids(gnomad_union)

        # Convert to pandas, process, and write outputs
        print(f"Processing variants for {gene_name}", flush=True)
        df = gnomad_union.to_pandas(flatten=True)
        df = process_dataframe(df, gene_name)
        write_outputs(df, gene_name)
        print(f"Successfully wrote outputs for {gene_name}", flush=True)

    except Exception as e:
        print(f"ERROR: {str(e)}", flush=True)
        import traceback
        traceback.print_exc()
        # Write empty output file to prevent Cromwell failure
        with open("gnomadcount.txt", 'w') as f:
            f.write("0")
        raise


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: get_gnomad_variants.py <gene_name> <chr_id> <start_locus> <end_locus> <hail_mem_gb>")
        sys.exit(1)

    gene_name = sys.argv[1]
    chr_id = sys.argv[2]
    start_locus = int(sys.argv[3])
    end_locus = int(sys.argv[4])
    hail_mem_gb = int(sys.argv[5])

    get_gnomad_variants(gene_name, chr_id, start_locus, end_locus, hail_mem_gb)
