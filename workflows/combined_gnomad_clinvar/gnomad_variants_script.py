import argparse
import os
import io
import json
import string
import pandas as pd
import re 

import gnomad
# Hail-specific packages
import hail as hl


parser = argparse.ArgumentParser(description='Get gnomad variants.')
parser.add_argument('-g', help='gene name', required=True)
parser.add_argument('-c', help='chromosome number', required=True)
parser.add_argument('-b', help='start position of gene', required=True)
parser.add_argument('-e', help='stop position of gene', required=True)
parser.add_argument('-o', default='gnomad_variants_MANE.csv',help='output CSV with gnomad variants')

args = parser.parse_args()

start_pos =hl.int32(args.b)
stop_pos =hl.int32(args.e)

from gnomad.resources.grch38.gnomad import public_release
v4exomes = public_release("exomes").ht()
v4exomes.count()


from gnomad.resources.grch38.gnomad import public_release
v4genomes = public_release("genomes").ht()
v4genomes.count()
#gnomad exome data
v4exomes_varid_sm = v4exomes.select( 
                                    v4exomes.freq,          #must be broken down further
                                    #v4exomes.joint_freq, 
                                    #v4exomes.vep.id,
                                    v4exomes.vep.allele_string, 
                                    v4exomes.vep.start, 
                                    v4exomes.vep.end, 
                                    #v4exomes.vep.strand, 
                                    v4exomes.vep.seq_region_name, 
                                    v4exomes.vep.variant_class,
                                    v4exomes.allele_info.allele_type,                                         
                                    v4exomes.in_silico_predictors.spliceai_ds_max,
                                    variant_type_exomes =   v4exomes.allele_info.variant_type,
                                    n_alt_alleles_exomes =   v4exomes.allele_info.n_alt_alleles,
                                    #txpt_allele_num = v4exomes.vep.transcript_consequences.allele_num,
                                    txpt_biotype = v4exomes.vep.transcript_consequences.biotype,
                                    txpt_consequence_terms = v4exomes.vep.transcript_consequences.consequence_terms,
                                    txpt_impact = v4exomes.vep.transcript_consequences.impact,
                                    #txpt_variant_allele = v4exomes.vep.transcript_consequences.variant_allele,
                                    
                                    txpt_amino_acids = v4exomes.vep.transcript_consequences.amino_acids,
                                    txpt_appris = v4exomes.vep.transcript_consequences.appris,
                                    txpt_canonical = v4exomes.vep.transcript_consequences.canonical,
                                    txpt_distance = v4exomes.vep.transcript_consequences.distance,
                                    txpt_domains = v4exomes.vep.transcript_consequences.domains,
                                    txpt_exon = v4exomes.vep.transcript_consequences.exon,
                                    txpt_hgvsc = v4exomes.vep.transcript_consequences.hgvsc,
                                    txpt_hgvsp = v4exomes.vep.transcript_consequences.hgvsp,
                                    txpt_hgvs_offset = v4exomes.vep.transcript_consequences.hgvs_offset,
                                    #txpt_gene_id = v4exomes.vep.transcript_consequences.gene_id,
                                    txpt_gene_pheno = v4exomes.vep.transcript_consequences.gene_pheno,
                                    txpt_gene_symbol = v4exomes.vep.transcript_consequences.gene_symbol,
                                    #txpt_gene_symbol_source = v4exomes.vep.transcript_consequences.gene_symbol_source,
                                    txpt_intron = v4exomes.vep.transcript_consequences.intron,
                                    txpt_lof = v4exomes.vep.transcript_consequences.lof,
                                    txpt_lof_flags = v4exomes.vep.transcript_consequences.lof_flags,
                                    txpt_lof_filter = v4exomes.vep.transcript_consequences.lof_filter,
                                    txpt_lof_info = v4exomes.vep.transcript_consequences.lof_info,
                                    txpt_mane_select = v4exomes.vep.transcript_consequences.mane_select,
                                    txpt_mane_plus_clinical = v4exomes.vep.transcript_consequences.mane_plus_clinical,
                                    txpt_protein_start = v4exomes.vep.transcript_consequences.protein_start,
                                    txpt_protein_end = v4exomes.vep.transcript_consequences.protein_end,
                                    #txpt_source = v4exomes.vep.transcript_consequences.source,
                                    txpt_transcript_id = v4exomes.vep.transcript_consequences.transcript_id,
                                    #txpt_tsl = v4exomes.vep.transcript_consequences.tsl,
                                    txpt_uniprot_isoform = v4exomes.vep.transcript_consequences.uniprot_isoform
                                

                                   )

#gnomad genome data
v4genomes_varid_sm = v4genomes.select(  #v4genomes.joint_freq,
                                      v4genomes.freq, 

                                                                        
                                    #v4genomes.vep.id,
                                    v4genomes.vep.allele_string, 
                                    v4genomes.vep.start, 
                                    v4genomes.vep.end, 
                                    #v4genomes.vep.strand, 
                                    v4genomes.vep.seq_region_name, 
                                    v4genomes.vep.variant_class,
                                    v4genomes.allele_info.variant_type,
                                    v4genomes.allele_info.allele_type,                                         
                                 
                                    v4genomes.in_silico_predictors.spliceai_ds_max,
                                    variant_type_genomes =   v4genomes.allele_info.variant_type,
                                    n_alt_alleles_genomes=   v4genomes.allele_info.n_alt_alleles,
                                    #txpt_allele_num = v4genomes.vep.transcript_consequences.allele_num,
                                    txpt_biotype = v4genomes.vep.transcript_consequences.biotype,
                                    txpt_consequence_terms = v4genomes.vep.transcript_consequences.consequence_terms,
                                    txpt_impact = v4genomes.vep.transcript_consequences.impact,
                                    #txpt_variant_allele = v4genomes.vep.transcript_consequences.variant_allele,
                                    
                                    txpt_amino_acids = v4genomes.vep.transcript_consequences.amino_acids,
                                    txpt_appris = v4genomes.vep.transcript_consequences.appris,
                                    txpt_canonical = v4genomes.vep.transcript_consequences.canonical,
                                    txpt_distance = v4genomes.vep.transcript_consequences.distance,
                                    txpt_domains = v4genomes.vep.transcript_consequences.domains,
                                    txpt_exon = v4genomes.vep.transcript_consequences.exon,
                                    txpt_hgvsc = v4genomes.vep.transcript_consequences.hgvsc,
                                    txpt_hgvsp = v4genomes.vep.transcript_consequences.hgvsp,
                                    txpt_hgvs_offset = v4genomes.vep.transcript_consequences.hgvs_offset,
                                    #txpt_gene_id = v4genomes.vep.transcript_consequences.gene_id,
                                    txpt_gene_pheno = v4genomes.vep.transcript_consequences.gene_pheno,
                                    txpt_gene_symbol = v4genomes.vep.transcript_consequences.gene_symbol,
                                    #txpt_gene_symbol_source = v4genomes.vep.transcript_consequences.gene_symbol_source,
                                    txpt_intron = v4genomes.vep.transcript_consequences.intron,
                                    txpt_lof = v4genomes.vep.transcript_consequences.lof,
                                    txpt_lof_flags = v4genomes.vep.transcript_consequences.lof_flags,
                                    txpt_lof_filter = v4genomes.vep.transcript_consequences.lof_filter,
                                    txpt_lof_info = v4genomes.vep.transcript_consequences.lof_info,
                                    txpt_mane_select = v4genomes.vep.transcript_consequences.mane_select,
                                    txpt_mane_plus_clinical = v4genomes.vep.transcript_consequences.mane_plus_clinical,
                                    txpt_protein_start = v4genomes.vep.transcript_consequences.protein_start,
                                    txpt_protein_end = v4genomes.vep.transcript_consequences.protein_end,
                                    #txpt_source = v4genomes.vep.transcript_consequences.source,
                                    txpt_transcript_id = v4genomes.vep.transcript_consequences.transcript_id,
                                    #txpt_tsl = v4genomes.vep.transcript_consequences.tsl,
                                    txpt_uniprot_isoform = v4genomes.vep.transcript_consequences.uniprot_isoform

                                        )

#The 'freq' annotation is an array, and each element of the array is a struct that contains the alternate allele count (AC), alternate allele frequency (AF), total number of alleles (AN), 
#and number of homozygous alternate individuals (homozygote_count) for a specific sample grouping.
#Use the 'freq_index_dict' global annotation to retrieve frequency information for a specific group of samples from the 'freq' array. 
#This global annotation is a dictionary keyed by sample grouping combinations whose values are the combination's index in the 'freq' array.


#freq, with ukb
v4exomes_varid_sm = v4exomes_varid_sm.annotate(exome_freq_main_adj_afr=v4exomes_varid_sm.freq[v4exomes_varid_sm.freq_index_dict['afr_adj']])

v4exomes_varid_sm = v4exomes_varid_sm.annotate(exome_freq_main_adj_amr=v4exomes_varid_sm.freq[v4exomes_varid_sm.freq_index_dict['amr_adj']])

v4exomes_varid_sm = v4exomes_varid_sm.annotate(exome_freq_main_adj_eas=v4exomes_varid_sm.freq[v4exomes_varid_sm.freq_index_dict['eas_adj']])

v4exomes_varid_sm = v4exomes_varid_sm.annotate(exome_freq_main_adj_nfe=v4exomes_varid_sm.freq[v4exomes_varid_sm.freq_index_dict['nfe_adj']])

v4exomes_varid_sm = v4exomes_varid_sm.annotate(exome_freq_main_adj_rmi=v4exomes_varid_sm.freq[v4exomes_varid_sm.freq_index_dict['remaining_adj']])

v4exomes_varid_sm = v4exomes_varid_sm.annotate(exome_freq_main_adj_sas=v4exomes_varid_sm.freq[v4exomes_varid_sm.freq_index_dict['sas_adj']])

v4exomes_varid_sm = v4exomes_varid_sm.annotate(exome_freq_main_adj_mid=v4exomes_varid_sm.freq[v4exomes_varid_sm.freq_index_dict['mid_adj']])

v4exomes_varid_sm = v4exomes_varid_sm.annotate(exome_freq_main_adj_fin=v4exomes_varid_sm.freq[v4exomes_varid_sm.freq_index_dict['fin_adj']])

v4exomes_varid_sm = v4exomes_varid_sm.annotate(exome_freq_main_adj_asj=v4exomes_varid_sm.freq[v4exomes_varid_sm.freq_index_dict['asj_adj']])



#freq, genomes, with ukb
v4genomes_varid_sm = v4genomes_varid_sm.annotate(genome_freq_adj_afr=v4genomes_varid_sm.freq[v4genomes_varid_sm.freq_index_dict['afr_adj']])

v4genomes_varid_sm = v4genomes_varid_sm.annotate(genome_freq_adj_amr=v4genomes_varid_sm.freq[v4genomes_varid_sm.freq_index_dict['amr_adj']])

v4genomes_varid_sm = v4genomes_varid_sm.annotate(genome_freq_adj_eas=v4genomes_varid_sm.freq[v4genomes_varid_sm.freq_index_dict['eas_adj']])

v4genomes_varid_sm = v4genomes_varid_sm.annotate(genome_freq_adj_nfe=v4genomes_varid_sm.freq[v4genomes_varid_sm.freq_index_dict['nfe_adj']])

v4genomes_varid_sm = v4genomes_varid_sm.annotate(genome_freq_adj_rmi=v4genomes_varid_sm.freq[v4genomes_varid_sm.freq_index_dict['remaining_adj']])

v4genomes_varid_sm = v4genomes_varid_sm.annotate(genome_freq_adj_sas=v4genomes_varid_sm.freq[v4genomes_varid_sm.freq_index_dict['sas_adj']])

v4genomes_varid_sm = v4genomes_varid_sm.annotate(genome_freq_adj_mid=v4genomes_varid_sm.freq[v4genomes_varid_sm.freq_index_dict['mid_adj']])

v4genomes_varid_sm = v4genomes_varid_sm.annotate(genome_freq_adj_fin=v4genomes_varid_sm.freq[v4genomes_varid_sm.freq_index_dict['fin_adj']])

v4genomes_varid_sm = v4genomes_varid_sm.annotate(genome_freq_adj_asj=v4genomes_varid_sm.freq[v4genomes_varid_sm.freq_index_dict['asj_adj']])


#The exomes and genomes databases do not contain gnomAD variant IDs. However, we can write a function that creates them from the locus and alleles data.
#the variant annotation functions
def get_ref_genome(locus: hl.expr.LocusExpression):
    """
    Expression for reference genome
    """
    #hl.genetics.Locus to get the Locus
    #hl.expr.LocusExpression to get the LocusExpression
    
    ref_gen = hl.str(gnomad.utils.reference_genome.get_reference_genome(locus).name)
    return ref_gen


def normalized_contig(contig: hl.expr.StringExpression) -> hl.expr.StringExpression:
    return hl.rbind(hl.str(contig).replace("^chr", ""), lambda c: hl.if_else(c == "MT", "M", c))

def chromosome_id(locus: hl.expr.LocusExpression, alleles: hl.expr.ArrayExpression, max_length: int = None):
    """
    Expression for computing chromosome of variant

    Args:
        max_length: (optional) length at which to truncate the <chrom> string

    Return:
        string: "<chrom>"
    """
    chr_id = hl.str(locus.contig) 

    if max_length is not None:
        return chr_id[0:max_length]

    return chr_id





def get_start_pos(locus: hl.expr.LocusExpression, alleles: hl.expr.ArrayExpression):
    """
    Expression for start position of both ref and alt allele
    """
    
    start_pos = hl.int32(locus.position)

    return start_pos

def get_end_pos_ref(locus: hl.expr.LocusExpression, alleles: hl.expr.ArrayExpression):
    """
    Expression for end position of ref  allele
    """
    ref_allele = hl.len(alleles[0]) 
    inmd = hl.int32(ref_allele -1)
    
    end_ref = hl.int32(locus.position) + inmd


    return end_ref

def get_end_pos_alt(locus: hl.expr.LocusExpression, alleles: hl.expr.ArrayExpression):
    """
    Expression for end position of ref  allele
    """
    
    alt_allele = hl.len(alleles[1])
    inmd = hl.int32(alt_allele -1)

    
    end_alt = hl.int32(locus.position) + inmd


    return end_alt


def variant_id_dash(locus: hl.expr.LocusExpression, alleles: hl.expr.ArrayExpression, max_length: int = None):
    """
    Expression for computing <chrom>-<pos>-<ref>-<alt>. Assumes alleles were split.

    Args:
        max_length: (optional) length at which to truncate the <chrom>-<pos>-<ref>-<alt> string

    Return:
        string: "<chrom>-<pos>-<ref>-<alt>"
    """
    refgen = get_ref_genome(locus)
    contig = normalized_contig(locus.contig)
    var_id = refgen + "-" + contig + "-" + hl.str(locus.position) + "-" + alleles[0] + "-" + alleles[1]

    if max_length is not None:
        return var_id[0:max_length]

    return var_id

def variant_id(locus: hl.expr.LocusExpression, alleles: hl.expr.ArrayExpression, max_length: int = None):
    """
    Expression for computing <chrom>:<pos>:<ref>:<alt>. Assumes alleles were split.

    Args:
        max_length: (optional) length at which to truncate the <chrom>:<pos>:<ref>:<alt> string

    Return:
        string: "<chrom>:<pos>:<ref>:<alt>"
    """
    refgen = get_ref_genome(locus)
    var_id = refgen + ":" + hl.str(locus.contig) + ":" + hl.str(locus.position) + ":" + alleles[0] + ":" + alleles[1]

    if max_length is not None:
        return var_id[0:max_length]

    return var_id

def variant_idCERFAC(locus: hl.expr.LocusExpression, alleles: hl.expr.ArrayExpression, max_length: int = None):
    """
    Expression for computing <chrom>:<pos>:<ref>:<alt>. Assumes alleles were split.

    Args:
        max_length: (optional) length at which to truncate the <chrom>:<pos>:<ref>:<alt> string

    Return:
        string: "<chrom>:<pos>:<ref>:<alt>"
    """
    refgen = get_ref_genome(locus)
    end_alt = get_end_pos_alt(locus,alleles)
    var_id = refgen + ":" + hl.str(locus.contig) + ":" + hl.str(locus.position) + ":" + hl.str(end_alt) + ":" + alleles[0] + ":" + alleles[1]

    if max_length is not None:
        return var_id[0:max_length]

    return var_id

def var_type_gnomad( alleles: hl.expr.ArrayExpression):
    """
    Expression for the allele type according to gnomad: 
    The possible return values are:
    "SNP"

    "MNP"

    "Insertion"

    "Deletion"

    "Complex"

    "Star"

    "Symbolic"

    "Unknown"


    """
    var_type = hl.str(hl.allele_type(hl.str(alleles[0]),hl.str(alleles[1])))

    return var_type




v4genomes_varid_sm = v4genomes_varid_sm.annotate(ref_genome=get_ref_genome(v4genomes_varid_sm.locus))


chr_string = "chr" + args.c
filtered_v4exomes = v4exomes_varid_sm.filter(v4exomes_varid_sm.locus.contig == chr_string )
filtered_v4genomes = v4genomes_varid_sm.filter(v4genomes_varid_sm.locus.contig == chr_string )
exomes_select = hl.filter_intervals(filtered_v4exomes, [hl.locus_interval(chr_string, start_pos, stop_pos, reference_genome='GRCh38')])
genomes_select = hl.filter_intervals(filtered_v4genomes, [hl.locus_interval(chr_string, start_pos, stop_pos, reference_genome='GRCh38')])

joined_gnomad_inner = genomes_select.join(exomes_select, how='inner')
joined_gnomad_inner = joined_gnomad_inner.annotate(exorgen="both")

joined_gnomad_inner = joined_gnomad_inner.drop(*(x for x in joined_gnomad_inner.row if re.search(r'_1', x)))

genomes_select_aj = genomes_select.anti_join(exomes_select)
genomes_select_aj = genomes_select_aj.annotate(exorgen="genomes")

exomes_select_aj = exomes_select.anti_join(genomes_select)
exomes_select_aj = exomes_select_aj.annotate(exorgen="exomes")


joined_gnomad_inner = joined_gnomad_inner.drop('freq')
exomes_select_aj = exomes_select_aj.drop('freq')
genomes_select_aj = genomes_select_aj.drop('freq')


#flattening allows for union later

joined_gnomad_inner = joined_gnomad_inner.flatten()
exomes_select_aj = exomes_select_aj.flatten()
genomes_select_aj = genomes_select_aj.flatten()

#adding the overlapping variants to the exomes exclusive variants
gnomad_union_1 = exomes_select_aj.union(joined_gnomad_inner, unify=True)
gnomad_union_1.count()

gnomad_union = gnomad_union_1.union(genomes_select_aj, unify=True)

#get reference genome column
gnomad_union = gnomad_union.annotate(ref_genome=get_ref_genome(gnomad_union.locus))

#chr id
gnomad_union = gnomad_union.annotate(chr_id=chromosome_id(gnomad_union.locus, gnomad_union.alleles))

#get start and end positions
#gnomad_union = gnomad_union.annotate(ref_start=get_start_pos(gnomad_union.locus, gnomad_union.alleles))
#gnomad_union = gnomad_union.annotate(ref_end=get_end_pos_ref(gnomad_union.locus, gnomad_union.alleles))

#gnomad_union = gnomad_union.annotate(alt_start=get_start_pos(gnomad_union.locus, gnomad_union.alleles))
#gnomad_union = gnomad_union.annotate(alt_end=get_end_pos_alt(gnomad_union.locus, gnomad_union.alleles))

#gnomad_union = gnomad_union.annotate(gnomad_variant_id=variant_id(gnomad_union.locus, gnomad_union.alleles))
#gnomad_union = gnomad_union.annotate(gnomad_variant_id_B=variant_id_dash(gnomad_union.locus, gnomad_union.alleles))

gnomad_union = gnomad_union.annotate(var_type_gnomad=var_type_gnomad(gnomad_union.alleles))
gnomad_union = gnomad_union.annotate(CERFAC_variant_id=variant_idCERFAC(gnomad_union.locus, gnomad_union.alleles))

gnomad_union_df = gnomad_union.to_pandas()
gnomad_union_df = gnomad_union_df.sort_index(axis=1)

badcols = [
       'txpt_amino_acids', 'txpt_appris', 'txpt_biotype', 'txpt_canonical',
       'txpt_consequence_terms', 'txpt_distance', 'txpt_domains', 'txpt_exon',
       'txpt_hgvsc','txpt_hgvsp','txpt_hgvs_offset',
       'txpt_gene_pheno', 'txpt_gene_symbol', 'txpt_impact', 'txpt_intron',
       'txpt_lof', 'txpt_lof_filter', 'txpt_lof_flags', 'txpt_lof_info',
       'txpt_mane_plus_clinical', 'txpt_mane_select', 'txpt_protein_end',
       'txpt_protein_start', 'txpt_transcript_id', 'txpt_uniprot_isoform']
gnomad_union_df = gnomad_union_df.explode(badcols)
gnomad_union_df = gnomad_union_df[gnomad_union_df.txpt_canonical == 1]
gnomad_union_df = gnomad_union_df[gnomad_union_df.txpt_gene_symbol == args.g]
gnomad_union_df = gnomad_union_df[gnomad_union_df['txpt_mane_select'].str.startswith('NM')]
#consequence terms  needs explanding again!
badcols = ['txpt_consequence_terms']
gnomad_union_df = gnomad_union_df.explode(badcols)
#should probably be optional along with splice variants
gnomad_union_df = gnomad_union_df[gnomad_union_df.txpt_consequence_terms != "upstream_gene_variant"]

#gnomad_union_df = gnomad_union_df[gnomad_union_df.txpt_consequence_terms != "intron_variant"]
#gnomad_union_df = gnomad_union_df[gnomad_union_df.txpt_consequence_terms != "splice_region_variant"]
#gnomad_union_df = gnomad_union_df[gnomad_union_df.txpt_consequence_terms != "splice_donor_region_variant"]
#gnomad_union_df = gnomad_union_df[gnomad_union_df.txpt_consequence_terms != "splice_donor_variant"]
#gnomad_union_df = gnomad_union_df[gnomad_union_df.txpt_consequence_terms != "splice_acceptor_variant"]
#gnomad_union_df = gnomad_union_df[gnomad_union_df.txpt_consequence_terms != "splice_donor_5th_base_variant"]
#gnomad_union_df = gnomad_union_df[gnomad_union_df.txpt_consequence_terms != "splice_polypyrimidine_tract_variant"]

#drop cols txpt_canonical, txpt_appris, txpt_biotype txpt_gene_pheno txpt_gene_symbol 'txpt_mane_plus_clinical', 'txpt_mane_select', seq region name, refgenome chrid
gnomad_union_df[['txpt_hgvsc' ]] = gnomad_union_df[['txpt_hgvsc' ]].astype('str')
gnomad_union_df['txpt_hgvsc'] = gnomad_union_df['txpt_hgvsc'].str.split(pat=":", n=1,  regex=False).str.get(1)

gnomad_union_df = gnomad_union_df.drop(columns=[ 'txpt_appris', 'txpt_biotype', 'txpt_canonical',
       'txpt_gene_pheno', 'txpt_gene_symbol', 
       'txpt_mane_plus_clinical', 'txpt_mane_select', 'txpt_protein_end',
       'txpt_protein_start', 'txpt_transcript_id'])

gnomad_union_df.to_csv( args.o,  sep=',', index=False )


