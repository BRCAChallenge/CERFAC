version 1.0

workflow get_gnomad_variants {

    meta {
        author: "Allison Cheney"
        email: "archeney@ucsc.edu"
        description: "Extract missense gnomad variants for a specified gene"
    }

    parameter_meta {
        GENE_NAME: "Name of relevant gene"
    }
    input {
        String GENE_NAME
        File clinvar_var
        
    }
    #The order in which the workflow block and task definitions are arranged in the script does not matter. 
    #Nor does the order of the call statements matter, as we'll see further on.
    call extract_gene_loc {
        input: GENE_NAME=GENE_NAME
    }
    call get_gnomad_variants{
        input: 
            GENE_NAME=GENE_NAME, 
            CHR_ID=extract_gene_loc.CHR_ID, 
            GENE_LENGTH=extract_gene_loc.GENE_LENGTH, 
            GENE_START_LOCUS=extract_gene_loc.GENE_START_LOCUS, 
            GENE_END_LOCUS=extract_gene_loc.GENE_END_LOCUS

    }
    call merge_variants{
        input:
            GENE_NAME=GENE_NAME,
            gnomadvar=get_gnomad_variants.gnomadvar,
            clinvar_var=clinvar_var
    }

    output{
        File output_calibration_variants_file = merge_variants.combined_var
        String output_gnomad_variants_count = get_gnomad_variants.gnomad_variants_count
        String output_total_variants_count = merge_variants.combined_variants_count
        Int gene_length = extract_gene_loc.GENE_LENGTH
    }

}




task extract_gene_loc {
    input {
        String GENE_NAME
        Int memSizeGB = 4
        Int threadCount = 1
        Int diskSizeGB = 5

    }

    command <<<
        set -eux -o pipefail

        esearch -db clinvar -query "~{GENE_NAME}[GENE] AND single_gene [PROP] AND homo sapiens [ORGN]" | efetch -format variationid -start 1 -stop 1 | 
        xtract -pattern VariationArchive  \
        -group ClassifiedRecord/SimpleAllele/GeneList/Gene/Location/SequenceLocation  -if SequenceLocation@Assembly -equals "GRCh38" -def "NA" \
            -element SequenceLocation@Assembly  SequenceLocation@Chr SequenceLocation@start  SequenceLocation@stop |
        tee gene_positions.txt

        awk -F '\t' 'NR == 1 {print $1}' gene_positions.txt | tee ASSEMBLY

        awk -F '\t' 'NR == 1 {print $2}' gene_positions.txt | tee CHR_ID
        
        if grep -q -m 1 "GRCh38" gene_positions.txt; then
            grep 'GRCh38' gene_positions.txt | awk -F '\t' '{print $3}' | tee GENE_START_LOCUS
            grep 'GRCh38' gene_positions.txt | awk -F '\t' '{print $4}' | tee GENE_END_LOCUS
        else 
            echo "couldn't find hg38 reference"
        fi
        pwd
    >>>

    output {
        File gene_results = "gene_positions.txt"
        String ASSEMBLY = read_string("ASSEMBLY")
        String CHR_ID = read_string("CHR_ID")
        Int GENE_START_LOCUS = read_int("GENE_START_LOCUS")
        Int GENE_END_LOCUS = read_int("GENE_END_LOCUS")
        Int GENE_LENGTH = GENE_END_LOCUS - GENE_START_LOCUS
        }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "allisoncheney/cerfac_terra:clinvar"
        preemptible: 1
    }

}


task get_gnomad_variants {
    input {
        Int memSizeGBbase = 15
        Int threadCount = 1
        Int diskSizeGBbase = 25
        String GENE_NAME
        String CHR_ID
        Int GENE_START_LOCUS
        Int GENE_END_LOCUS
        Int GENE_LENGTH
    }
    #private declaration here I guess
    Int overmilModifier = if GENE_LENGTH >= 1000000 then "30" else "0"
    Int overtwomilModifier = if GENE_LENGTH >= 2000000 then "45" else "0"

    Int memory_calc = memSizeGBbase + overmilModifier + overtwomilModifier
    Int hailMemSizeGB = floor((0.8* memory_calc )-2)
    Int diskSizeGB = diskSizeGBbase + overmilModifier + overtwomilModifier + 10


    command <<<
        set -ex -o pipefail
        echo "MEM_SIZE=$MEM_SIZE" >&2
        echo "MEM_UNIT=$MEM_UNIT" >&2


        python3 <<CODE

        import os
        import io
        import json
        import string
        import pandas as pd
        import re
        import logging

        import gnomad
        import hail as hl
        
        hl.init(spark_conf={'spark.driver.memory': "~{hailMemSizeGB}g"})


        
        start_pos =hl.int32("~{GENE_START_LOCUS}")
        stop_pos =hl.int32("~{GENE_END_LOCUS}")

        from gnomad.resources.grch38.gnomad import public_release
        v4exomes = public_release("exomes").ht()
        v4exomes.count()


        from gnomad.resources.grch38.gnomad import public_release
        v4genomes = public_release("genomes").ht()
        v4genomes.count()




        v4exomes_varid_sm = v4exomes.select(v4exomes.freq, 
                                            v4exomes.vep.allele_string, 
                                            v4exomes.vep.start, 
                                            v4exomes.vep.end, 
                                            v4exomes.vep.seq_region_name, 
                                            v4exomes.vep.variant_class,
                                            v4exomes.in_silico_predictors.spliceai_ds_max,
                                            n_alt_alleles_exomes =   v4exomes.allele_info.n_alt_alleles,
                                            txpt_biotype = v4exomes.vep.transcript_consequences.biotype,
                                            variant_effect = v4exomes.vep.transcript_consequences.consequence_terms,
                                            txpt_impact = v4exomes.vep.transcript_consequences.impact,
                                            VRS_Allele_IDs = v4exomes.info.vrs.VRS_Allele_IDs,
                                            allele_list_VRS = v4exomes.info.vrs.VRS_States,
                                            pos_VRS_starts = v4exomes.info.vrs.VRS_Starts,
                                            pos_VRS_stops = v4exomes.info.vrs.VRS_Ends,

                                            txpt_amino_acids = v4exomes.vep.transcript_consequences.amino_acids,
                                            txpt_appris = v4exomes.vep.transcript_consequences.appris,
                                            txpt_canonical = v4exomes.vep.transcript_consequences.canonical,
                                            txpt_distance = v4exomes.vep.transcript_consequences.distance,
                                            txpt_domains = v4exomes.vep.transcript_consequences.domains,
                                            txpt_exon = v4exomes.vep.transcript_consequences.exon,
                                            txpt_hgvsc = v4exomes.vep.transcript_consequences.hgvsc,
                                            txpt_hgvsp = v4exomes.vep.transcript_consequences.hgvsp,
                                            txpt_gene_pheno = v4exomes.vep.transcript_consequences.gene_pheno,
                                            txpt_gene_symbol = v4exomes.vep.transcript_consequences.gene_symbol,
                                            txpt_intron = v4exomes.vep.transcript_consequences.intron,
                                            txpt_lof = v4exomes.vep.transcript_consequences.lof,
                                            txpt_lof_flags = v4exomes.vep.transcript_consequences.lof_flags,
                                            txpt_lof_filter = v4exomes.vep.transcript_consequences.lof_filter,
                                            txpt_lof_info = v4exomes.vep.transcript_consequences.lof_info,
                                            txpt_mane_select = v4exomes.vep.transcript_consequences.mane_select,
                                            txpt_mane_plus_clinical = v4exomes.vep.transcript_consequences.mane_plus_clinical,
                                            txpt_protein_start = v4exomes.vep.transcript_consequences.protein_start,
                                            txpt_protein_end = v4exomes.vep.transcript_consequences.protein_end,
                                            txpt_transcript_id = v4exomes.vep.transcript_consequences.transcript_id,
                                            txpt_uniprot_isoform = v4exomes.vep.transcript_consequences.uniprot_isoform 
                                            )
        v4genomes_varid_sm = v4genomes.select( v4genomes.freq, 
                                            v4genomes.vep.allele_string, 
                                            v4genomes.vep.start, 
                                            v4genomes.vep.end, 
                                            v4genomes.vep.seq_region_name, 
                                            v4genomes.vep.variant_class,
                                            v4genomes.in_silico_predictors.spliceai_ds_max,
                                            VRS_Allele_IDs = v4genomes.info.vrs.VRS_Allele_IDs,
                                            allele_list_VRS = v4genomes.info.vrs.VRS_States,
                                            pos_VRS_starts = v4genomes.info.vrs.VRS_Starts,
                                            pos_VRS_stops = v4genomes.info.vrs.VRS_Ends,
                                            n_alt_alleles_genomes=   v4genomes.allele_info.n_alt_alleles,
                                            txpt_biotype = v4genomes.vep.transcript_consequences.biotype,
                                            variant_effect = v4genomes.vep.transcript_consequences.consequence_terms,
                                            txpt_impact = v4genomes.vep.transcript_consequences.impact,
                                            
                                            txpt_amino_acids = v4genomes.vep.transcript_consequences.amino_acids,
                                            txpt_appris = v4genomes.vep.transcript_consequences.appris,
                                            txpt_canonical = v4genomes.vep.transcript_consequences.canonical,
                                            txpt_distance = v4genomes.vep.transcript_consequences.distance,
                                            txpt_domains = v4genomes.vep.transcript_consequences.domains,
                                            txpt_exon = v4genomes.vep.transcript_consequences.exon,
                                            txpt_hgvsc = v4genomes.vep.transcript_consequences.hgvsc,
                                            txpt_hgvsp = v4genomes.vep.transcript_consequences.hgvsp,
                                            txpt_gene_pheno = v4genomes.vep.transcript_consequences.gene_pheno,
                                            txpt_gene_symbol = v4genomes.vep.transcript_consequences.gene_symbol,
                                            txpt_intron = v4genomes.vep.transcript_consequences.intron,
                                            txpt_lof = v4genomes.vep.transcript_consequences.lof,
                                            txpt_lof_flags = v4genomes.vep.transcript_consequences.lof_flags,
                                            txpt_lof_filter = v4genomes.vep.transcript_consequences.lof_filter,
                                            txpt_lof_info = v4genomes.vep.transcript_consequences.lof_info,
                                            txpt_mane_select = v4genomes.vep.transcript_consequences.mane_select,
                                            txpt_mane_plus_clinical = v4genomes.vep.transcript_consequences.mane_plus_clinical,
                                            txpt_protein_start = v4genomes.vep.transcript_consequences.protein_start,
                                            txpt_protein_end = v4genomes.vep.transcript_consequences.protein_end,
                                            txpt_transcript_id = v4genomes.vep.transcript_consequences.transcript_id,
                                            txpt_uniprot_isoform = v4genomes.vep.transcript_consequences.uniprot_isoform

                                                )

        v4exomes_varid_sm = v4exomes_varid_sm.annotate(exome_freq_main_adj_afr=v4exomes_varid_sm.freq[v4exomes_varid_sm.freq_index_dict['afr_adj']])

        v4exomes_varid_sm = v4exomes_varid_sm.annotate(exome_freq_main_adj_amr=v4exomes_varid_sm.freq[v4exomes_varid_sm.freq_index_dict['amr_adj']])

        v4exomes_varid_sm = v4exomes_varid_sm.annotate(exome_freq_main_adj_eas=v4exomes_varid_sm.freq[v4exomes_varid_sm.freq_index_dict['eas_adj']])

        v4exomes_varid_sm = v4exomes_varid_sm.annotate(exome_freq_main_adj_nfe=v4exomes_varid_sm.freq[v4exomes_varid_sm.freq_index_dict['nfe_adj']])

        v4exomes_varid_sm = v4exomes_varid_sm.annotate(exome_freq_main_adj_rmi=v4exomes_varid_sm.freq[v4exomes_varid_sm.freq_index_dict['remaining_adj']])

        v4exomes_varid_sm = v4exomes_varid_sm.annotate(exome_freq_main_adj_sas=v4exomes_varid_sm.freq[v4exomes_varid_sm.freq_index_dict['sas_adj']])

        v4exomes_varid_sm = v4exomes_varid_sm.annotate(exome_freq_main_adj_mid=v4exomes_varid_sm.freq[v4exomes_varid_sm.freq_index_dict['mid_adj']])

        v4exomes_varid_sm = v4exomes_varid_sm.annotate(exome_freq_main_adj_fin=v4exomes_varid_sm.freq[v4exomes_varid_sm.freq_index_dict['fin_adj']])

        v4exomes_varid_sm = v4exomes_varid_sm.annotate(exome_freq_main_adj_asj=v4exomes_varid_sm.freq[v4exomes_varid_sm.freq_index_dict['asj_adj']])



        v4genomes_varid_sm = v4genomes_varid_sm.annotate(genome_freq_adj_afr=v4genomes_varid_sm.freq[v4genomes_varid_sm.freq_index_dict['afr_adj']])

        v4genomes_varid_sm = v4genomes_varid_sm.annotate(genome_freq_adj_amr=v4genomes_varid_sm.freq[v4genomes_varid_sm.freq_index_dict['amr_adj']])

        v4genomes_varid_sm = v4genomes_varid_sm.annotate(genome_freq_adj_eas=v4genomes_varid_sm.freq[v4genomes_varid_sm.freq_index_dict['eas_adj']])

        v4genomes_varid_sm = v4genomes_varid_sm.annotate(genome_freq_adj_nfe=v4genomes_varid_sm.freq[v4genomes_varid_sm.freq_index_dict['nfe_adj']])

        v4genomes_varid_sm = v4genomes_varid_sm.annotate(genome_freq_adj_rmi=v4genomes_varid_sm.freq[v4genomes_varid_sm.freq_index_dict['remaining_adj']])

        v4genomes_varid_sm = v4genomes_varid_sm.annotate(genome_freq_adj_sas=v4genomes_varid_sm.freq[v4genomes_varid_sm.freq_index_dict['sas_adj']])

        v4genomes_varid_sm = v4genomes_varid_sm.annotate(genome_freq_adj_mid=v4genomes_varid_sm.freq[v4genomes_varid_sm.freq_index_dict['mid_adj']])

        v4genomes_varid_sm = v4genomes_varid_sm.annotate(genome_freq_adj_fin=v4genomes_varid_sm.freq[v4genomes_varid_sm.freq_index_dict['fin_adj']])

        v4genomes_varid_sm = v4genomes_varid_sm.annotate(genome_freq_adj_asj=v4genomes_varid_sm.freq[v4genomes_varid_sm.freq_index_dict['asj_adj']])


        def get_ref_genome(locus: hl.expr.LocusExpression):
            """
            Expression for reference genome
            """
            
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

        def get_alt_allele_len(locus: hl.expr.LocusExpression, alleles: hl.expr.ArrayExpression):
            """
            Expression alt allele
            """
            
            alt_allele = hl.len(alleles[1])



            return alt_allele
        def get_ref_allele_len(locus: hl.expr.LocusExpression, alleles: hl.expr.ArrayExpression):
            """
            Expression alt allele
            """
            
            ref_allele = hl.len(alleles[0])



            return ref_allele
        def get_alt_allele(locus: hl.expr.LocusExpression, alleles: hl.expr.ArrayExpression):
            """
            Expression alt allele
            """
            
            alt_allele = hl.str(alleles[1])



            return alt_allele

        def get_ref_allele(locus: hl.expr.LocusExpression, alleles: hl.expr.ArrayExpression):
            """
            Expression alt allele
            """
            
            ref_allele = hl.str(alleles[0])



            return ref_allele
        def get_VRS_start_ref(pos_VRS_starts: hl.expr.ArrayNumericExpression):
            """
            Position ref allele
            """
            
            ref_pos_start = hl.str(pos_VRS_starts[0])



            return ref_pos_start
        def get_VRS_start_alt(locus: hl.expr.LocusExpression, pos_VRS_starts: hl.expr.ArrayNumericExpression):
            """
            Position alt allele
            """
            
            alt_pos_start = hl.str(pos_VRS_starts[1])



            return alt_pos_start
        def get_VRS_stop_ref(pos_VRS_stops: hl.expr.ArrayNumericExpression):
            """
            Position ref allele
            """
            
            ref_pos_stop = hl.str(pos_VRS_stops[0])



            return ref_pos_stop
        def get_VRS_stop_alt(pos_VRS_stops: hl.expr.ArrayNumericExpression):
            """
            Position alt allele
            """
            
            alt_pos_stop = hl.str(pos_VRS_stops[1])



            return alt_pos_stop

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

        def variant_idCERFAC_orig(locus: hl.expr.LocusExpression, alleles: hl.expr.ArrayExpression, max_length: int = None):
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

        def variant_idCERFAC_VCF(locus: hl.expr.LocusExpression, alleles: hl.expr.ArrayExpression, max_length: int = None):
            """
            Expression for computing <chrom>:<pos>:<ref>:<alt>. Assumes alleles were split.

            Args:
                max_length: (optional) length at which to truncate the <chrom>:<pos>:<ref>:<alt> string

            Return:
                string: "<chrom>-<pos>-<ref>-<alt>"
            """
            refgen = get_ref_genome(locus)
            end_alt = get_end_pos_alt(locus,alleles)
            var_id =  hl.str(locus.contig) + "-" + hl.str(locus.position) + "-" + alleles[0] + "-" + alleles[1]

            if max_length is not None:
                return var_id[0:max_length]

            return var_id




        v4genomes_varid_sm = v4genomes_varid_sm.annotate(ref_genome=get_ref_genome(v4genomes_varid_sm.locus))


        chr_string = "chr" + "~{CHR_ID}"
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



        joined_gnomad_inner = joined_gnomad_inner.flatten()
        exomes_select_aj = exomes_select_aj.flatten()
        genomes_select_aj = genomes_select_aj.flatten()

        gnomad_union_1 = exomes_select_aj.union(joined_gnomad_inner, unify=True)
        gnomad_union_1.count()

        gnomad_union = gnomad_union_1.union(genomes_select_aj, unify=True)

        gnomad_union = gnomad_union.annotate(ref_genome=get_ref_genome(gnomad_union.locus))

        gnomad_union = gnomad_union.annotate(chr_id=chromosome_id(gnomad_union.locus, gnomad_union.alleles))


        gnomad_union = gnomad_union.annotate(CERFAC_variant_id_VCF=variant_idCERFAC_VCF(gnomad_union.locus, gnomad_union.alleles))
        gnomad_union = gnomad_union.annotate(allele_ref=get_ref_allele(gnomad_union.locus, gnomad_union.alleles))
        gnomad_union = gnomad_union.annotate(allele_alt=get_alt_allele(gnomad_union.locus, gnomad_union.alleles))
        gnomad_union = gnomad_union.annotate(variant_length_ref=get_ref_allele_len(gnomad_union.locus, gnomad_union.alleles))
        gnomad_union = gnomad_union.annotate(variant_length_alt=get_alt_allele_len(gnomad_union.locus, gnomad_union.alleles))

        gnomad_union = gnomad_union.annotate(pos_start_alt_vrs=get_VRS_start_alt(gnomad_union.locus, gnomad_union.pos_VRS_starts))
        gnomad_union = gnomad_union.annotate(pos_start_ref_vrs=get_VRS_start_ref(gnomad_union.pos_VRS_starts))

        gnomad_union = gnomad_union.annotate(pos_stop_alt_vrs=get_VRS_stop_alt(gnomad_union.pos_VRS_stops))
        gnomad_union = gnomad_union.annotate(pos_stop_ref_vrs=get_VRS_stop_ref(gnomad_union.pos_VRS_stops))




        gnomad_union_df = gnomad_union.to_pandas(flatten=True)
        gnomad_union_df = gnomad_union_df.sort_index(axis=1)

        badcols = [
            'txpt_amino_acids', 'txpt_appris', 'txpt_biotype', 'txpt_canonical',
            'variant_effect', 'txpt_distance', 'txpt_domains', 'txpt_exon',
            'txpt_hgvsc','txpt_hgvsp',
            'txpt_gene_pheno', 'txpt_gene_symbol', 'txpt_impact', 'txpt_intron',
            'txpt_lof', 'txpt_lof_filter', 'txpt_lof_flags', 'txpt_lof_info',
            'txpt_mane_plus_clinical', 'txpt_mane_select', 'txpt_protein_end',
            'txpt_protein_start', 'txpt_transcript_id', 'txpt_uniprot_isoform']
        gnomad_union_df = gnomad_union_df.explode(badcols)
        gnomad_union_df = gnomad_union_df[gnomad_union_df.txpt_canonical == 1]
        gnomad_union_df = gnomad_union_df.dropna(subset=['txpt_gene_symbol'])
        gnomad_union_df = gnomad_union_df[gnomad_union_df.txpt_gene_symbol == "~{GENE_NAME}"]
        gnomad_union_df = gnomad_union_df[gnomad_union_df['txpt_mane_select'].str.startswith('NM')]

        gnomad_union_df[['variant_effect' ]] = gnomad_union_df[['variant_effect' ]].astype('str')
        gnomad_union_df.loc[:, ('variant_effect')] = gnomad_union_df.loc[:, ('variant_effect')].astype('str').str.replace('[','')
        gnomad_union_df.loc[:, ('variant_effect')] = gnomad_union_df.loc[:, ('variant_effect')].astype('str').str.replace(']','')
        gnomad_union_df.loc[:, ('variant_effect')] = gnomad_union_df.loc[:, ('variant_effect')].astype('str').str.replace("'","")



        badcols = ['txpt_uniprot_isoform']
        gnomad_union_df = gnomad_union_df.explode(badcols)
        gnomad_union_df = gnomad_union_df[gnomad_union_df.variant_effect != "upstream_gene_variant"]
        gnomad_union_df = gnomad_union_df[gnomad_union_df.variant_effect != "5_prime_UTR_variant"]
        gnomad_union_df = gnomad_union_df[gnomad_union_df.variant_effect != "3_prime_UTR_variant"]
        gnomad_union_df = gnomad_union_df[gnomad_union_df.variant_effect != "intron_variant"]


        gnomad_union_df[['txpt_hgvsc' ]] = gnomad_union_df[['txpt_hgvsc' ]].astype('str')
        gnomad_union_df['CERFAC_variant_id_HGVS_long'] = gnomad_union_df[['ref_genome', 'locus','txpt_hgvsc' ]].astype(str).agg(':'.join, axis=1)

        gnomad_union_df['txpt_hgvsc_short'] = gnomad_union_df['txpt_hgvsc'].str.split(pat=":", n=1,  regex=False).str.get(1)
        gnomad_union_df['hgvs_pro'] = gnomad_union_df['txpt_hgvsp'].str.split(pat=":", n=1,  regex=False).str.get(1)

        gnomad_union_df['CERFAC_variant_id_HGVS_short'] = gnomad_union_df[['ref_genome','locus','txpt_hgvsc_short' ]].astype(str).agg(':'.join, axis=1)
        gnomad_union_df['CERFAC_variant_id_VCF'] = gnomad_union_df[['CERFAC_variant_id_VCF' ]].astype(str).replace("chr", "")


        gnomad_union_df = gnomad_union_df.drop(columns=[ 'txpt_appris',  'txpt_distance', 'txpt_biotype', 'txpt_canonical',
            'txpt_gene_pheno', 'txpt_gene_symbol', 'seq_region_name','chr_id','ref_genome',
            'txpt_protein_end', 'txpt_protein_start'])

        gnomad_union_df[['txpt_exon' ]] = gnomad_union_df[['txpt_exon' ]].astype('str')
        gnomad_union_df[['txpt_intron' ]] = gnomad_union_df[['txpt_intron' ]].astype('str')

        gnomad_variants_count_pd = str(gnomad_union_df['txpt_hgvsc_short'].nunique())
        file_name = "gnomadcount.txt"
        with open(file_name, 'w') as x_file:
            x_file.write(gnomad_variants_count_pd)


        gnomad_union_df = gnomad_union_df.rename(columns={"alleles": "allele_list",   "locus": "pos_VCF",    "exorgen": "set"}, errors='raise')
        gnomad_union_df = gnomad_union_df.rename(columns={"start": "pos_start_vep",   "end": "pos_stop_vep"}, errors='raise')
        gnomad_union_df = gnomad_union_df.sort_index(axis=1)
        gnomad_union_df['variant_source']="gnomAD"
        gnomad_union_df = gnomad_union_df.add_suffix('_gnomad')
        gnomad_union_df['hgvs_cdna_gnomad'] = gnomad_union_df[['txpt_mane_select_gnomad', 'txpt_hgvsc_short_gnomad' ]].astype(str).agg(':'.join, axis=1)

        


        


        gnomad_union_df.to_csv( "gnomad_variants_MANE.csv",  sep=',', index=False )
        CODE


    >>>

    output {
        File gnomadvar = "gnomad_variants_MANE.csv"
        String gnomad_variants_count = read_string("gnomadcount.txt")
    }

    runtime {
        memory: memory_calc + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "allisoncheney/cerfac_terra:gnomad"
        maxRetries: 0
        preemptible: 1
    }
}



task merge_variants {
    input {
        File gnomadvar
        File clinvar_var
        Int memSizeGB = 6
        Int threadCount = 1
        Int diskSizeGB = 5*round(size(gnomadvar, "GB")) + 10
        String GENE_NAME
    }

    command <<<
        set -eux -o pipefail
        python3 <<CODE
        import os
        import io
        import json
        import string
        import pandas as pd
        import numpy as np



        cv_table = pd.read_csv("~{clinvar_var}", sep=',' )
        gnomad_vars = pd.read_csv("~{gnomadvar}", sep=',' )
        gnomad_vars = gnomad_vars.rename(columns={"txpt_hgvsc_short_gnomad": "hgvs_nt"}, errors='raise')
        cv_table = cv_table.rename(columns={"txpt_hgvsc_clinvar": "hgvs_nt"}, errors='raise')


        combined = gnomad_vars.set_index('hgvs_nt').join(cv_table.set_index('hgvs_nt'), how='outer', lsuffix='_gnomad', rsuffix='_clinvar' )


        combined.sort_values(['hgvs_nt'])
        combined = combined.reset_index()


        combined_variants_count_pd = str(combined['hgvs_nt'].nunique())


        file_name = "combinedcount.txt"
        with open(file_name, 'w') as x_file:
            x_file.write(combined_variants_count_pd)
        combined = combined.set_index('hgvs_nt')
        combined = combined.sort_index(axis=1)
        combined['variant_source'] = combined.variant_source_gnomad.astype(str).str.cat(combined.variant_source_clinvar.astype(str), sep='', na_rep=None)
        combined['variant_source'] = combined['variant_source'].replace(to_replace=dict(gnomADClinVar="gnomAD and ClinVar", nanClinVar="ClinVar only", gnomADnan="gnomAD only"))

        rearrcols = [  'txpt_hgvsc_from_ID_clinvar', 'hgvs_pro_clinvar', 'hgvs_pro_gnomad',  'txpt_hgvsp_gnomad', 'variant_source',  'set_gnomad',   'n_alt_alleles_exomes_gnomad',   'n_alt_alleles_genomes_gnomad',   
        'number_submissions_clinvar', 
        'overall_germline_classification_clinvar','submission_germline_classification_clinvar',
        'overall_oncogenicity_classification_clinvar','submission_oncogenicity_classification_clinvar',
        'overall_somatic_classification_clinvar','submission_somatic_classification_clinvar',
        'functional_category_clinvar','functional_comment_clinvar','functional_result_clinvar',
        'variant_class_gnomad','variant_class_clinvar', 'variant_effect_gnomad','variant_effect_clinvar',
        'MG_disease_name_clinvar','ContributesToAggregateClassification_clinvar',
        'comment_clinvar',
        'VCV_ID_clinvar',
        'ClinVar_variant_ID_clinvar', 'SCV_ID_clinvar']

        combined = combined[rearrcols + [c for c in combined.columns if c not in rearrcols]]
        combined = combined.drop(columns=[ 'variant_source_clinvar',  'variant_source_gnomad' ])
        combined = combined.reset_index()
        from natsort import index_natsorted
        combined = combined.sort_values( by='hgvs_nt', key=lambda x: np.argsort(index_natsorted(combined["hgvs_nt"])))

        combined = combined.set_index('hgvs_nt')







        combined.to_csv( "~{GENE_NAME}_calibration_variants.csv",  sep=',', index=True )



        CODE

    >>>

    output {
        File combined_var = "~{GENE_NAME}_calibration_variants.csv"
        String combined_variants_count = read_string("combinedcount.txt")
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "allisoncheney/cerfac_terra:clinvar"
        maxRetries: 3
        preemptible: 1
    }
}





