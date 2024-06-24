version 1.0

workflow annotate_functional_variants {

    meta {
        author: "Allison Cheney"
        email: "archeney@ucsc.edu"
        description: "Extract variants for a specified gene"
    }

    parameter_meta {
        GENE_NAME: "Name of relevant gene"
    }
    input {
        String GENE_NAME
        
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
            GENE_START_LOCUS=extract_gene_loc.GENE_START_LOCUS, 
            GENE_END_LOCUS=extract_gene_loc.GENE_END_LOCUS

    }
    call  get_clinvar_variants_file{
        input: GENE_NAME=GENE_NAME
    }
    call extract_clinvar_variants_traitmap{
        input: 
            GENE_NAME=GENE_NAME,
            basicxml=get_clinvar_variants_file.basicxml
    }
    call extract_clinvar_variants_traitset{
        input: 
            GENE_NAME=GENE_NAME,
            basicxml=get_clinvar_variants_file.basicxml
    }
    call extract_clinvar_variants_basic{
        input: 
            GENE_NAME=GENE_NAME,
            basicxml=get_clinvar_variants_file.basicxml
    }
    call merge_clinvar_variants{
        input: 
            basiccv=extract_clinvar_variants_basic.basiccv, 
            traitset=extract_clinvar_variants_traitset.traitset, 
            traitmap=extract_clinvar_variants_traitmap.traitmap
    }
    call merge_variants{
        input:
            GENE_NAME=GENE_NAME,
            gnomadvar=get_gnomad_variants.gnomadvar,
            clinvar_var=merge_clinvar_variants.clinvar_var
    }


    output{
        File combined_var = merge_variants.combined_var
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
        Int memSizeGB = 10
        Int threadCount = 1
        Int diskSizeGB = 20
        String GENE_NAME
        String CHR_ID
        Int GENE_START_LOCUS
        Int GENE_END_LOCUS
    }

    command <<<
        set -eux -o pipefail

        python3 <<CODE

        import os
        import io
        import json
        import string
        import pandas as pd
        import re 

        import gnomad
        # Hail-specific packages
        import hail as hl

        start_pos =hl.int32("~{GENE_START_LOCUS}")
        stop_pos =hl.int32("~{GENE_END_LOCUS}")

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
                                            VRS_starts = v4exomes.info.vrs.VRS_Starts,
                                            VRS_stops = v4exomes.info.vrs.VRS_Ends,

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

                                            VRS_starts = v4genomes.info.vrs.VRS_Starts,
                                            VRS_stops = v4genomes.info.vrs.VRS_Ends,
                                            
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
                string: "<chrom>:<pos>:<ref>:<alt>"
            """
            refgen = get_ref_genome(locus)
            end_alt = get_end_pos_alt(locus,alleles)
            var_id = refgen + ":" + hl.str(locus.contig) + ":" + hl.str(locus.position) + ":" + alleles[0] + ":" + alleles[1]

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
        gnomad_union = gnomad_union.annotate(CERFAC_variant_id_orig=variant_idCERFAC_orig(gnomad_union.locus, gnomad_union.alleles))

        gnomad_union = gnomad_union.annotate(CERFAC_variant_id_VCF=variant_idCERFAC_VCF(gnomad_union.locus, gnomad_union.alleles))

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
        gnomad_union_df = gnomad_union_df[gnomad_union_df.txpt_gene_symbol == "~{GENE_NAME}"]
        #there are duplicate rows for this column with either the ENS ID or the NM ID, choose the NM ID
        gnomad_union_df = gnomad_union_df[gnomad_union_df['txpt_mane_select'].str.startswith('NM')]


        #consequence terms  needs explanding again!
        badcols = ['txpt_consequence_terms']
        gnomad_union_df = gnomad_union_df.explode(badcols)

        #commenting out because there's more than one start and stop causing duplicate rows in gnomad
        #badcols = [ 'VRS_starts','VRS_stops']
        #gnomad_union_df = gnomad_union_df.explode(badcols)

        badcols = ['txpt_uniprot_isoform']
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

        gnomad_union_df[['txpt_hgvsc' ]] = gnomad_union_df[['txpt_hgvsc' ]].astype('str')
        gnomad_union_df['CERFAC_variant_id_HGVS_long'] = gnomad_union_df[['ref_genome', 'locus','txpt_hgvsc' ]].astype(str).agg(':'.join, axis=1)

        gnomad_union_df['txpt_hgvsc'] = gnomad_union_df['txpt_hgvsc'].str.split(pat=":", n=1,  regex=False).str.get(1)

        gnomad_union_df['CERFAC_variant_id_HGVS_short'] = gnomad_union_df[['ref_genome','locus','txpt_hgvsc' ]].astype(str).agg(':'.join, axis=1)


        #drop cols txpt_canonical, txpt_appris, txpt_biotype txpt_gene_pheno txpt_gene_symbol 'txpt_mane_plus_clinical', 'txpt_mane_select', seq region name, refgenome chrid


        #probably remove 'variant_type', 'variant_type_exomes','variant_type_genomes'
        gnomad_union_df = gnomad_union_df.drop(columns=[ 'txpt_appris',  'txpt_distance','txpt_hgvs_offset', 'txpt_biotype', 'txpt_canonical',
            'txpt_gene_pheno', 'txpt_gene_symbol', 'seq_region_name','chr_id','ref_genome',
            # 'txpt_mane_plus_clinical', 'txpt_mane_select',  'txpt_transcript_id',
            'txpt_protein_end', 'txpt_protein_start'])
        gnomad_union_df.to_csv( "gnomad_variants_MANE.csv",  sep=',', index=False )
        CODE


    >>>

    output {
        File gnomadvar = "gnomad_variants_MANE.csv"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "allisoncheney/cerfac_terra:gnomad"
        maxRetries: 3
        preemptible: 1
    }
}


task get_clinvar_variants_file {
    input {
        String GENE_NAME
        Int memSizeGB = 10
        Int threadCount = 1
        Int diskSizeGB = 10
    }

    command <<<
        set -eux -o pipefail

        esearch -db clinvar -query "~{GENE_NAME}[GENE] AND single_gene [PROP] AND homo sapiens [ORGN] AND (varlen 49 or less[FILTER]) NOT (near gene upstream[PROP]) NOT (near gene downstream[PROP])" |
        efetch -format variationid > basic.xml
    >>>

    output {
        File basicxml  = "basic.xml"
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

task extract_clinvar_variants_traitmap {
    input {
        String GENE_NAME
        File basicxml  
        Int memSizeGB = 10
        Int threadCount = 1
        Int diskSizeGB = 5*round(size(basicxml, "GB")) + 2

    }

    command <<<
        set -eux -o pipefail

        cat  ~{basicxml} |
        xtract -pattern VariationArchive -def 'NA' -KEYVCV VariationArchive@Accession \
                -group TraitMapping -deq '\n' -def 'None' -lbl 'traitmapping' -element '&KEYVCV' @ClinicalAssertionID @TraitType MedGen@CUI MedGen@Name > ~{GENE_NAME}_traitmapping.txt

    >>>

    output {
        File traitmap  = "~{GENE_NAME}_traitmapping.txt"
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


task extract_clinvar_variants_traitset {
    input {
        String GENE_NAME
        File basicxml  
        Int memSizeGB = 10
        Int threadCount = 1
        Int diskSizeGB = 5*round(size(basicxml, "GB")) + 2

    }

    command <<<
        set -eux -o pipefail

        cat  ~{basicxml} |
        xtract -pattern VariationArchive -def 'NA' -KEYVCV VariationArchive@Accession \
            -group GermlineClassification/ConditionList \
                -block TraitSet -deq '\n' -def 'NA'   -TSID TraitSet@ID  -CONTRIB TraitSet@ContributesToAggregateClassification  -TSTYPE  TraitSet@Type \
                    -section Trait -deq '\n' -def 'NA'  -element '&KEYVCV'  '&TSID' '&TSTYPE' Trait@ID Trait@Type  '&CONTRIB'  \
                        -subset Trait/XRef -if XRef@DB -equals 'MedGen'  -element   XRef@ID  \
            -group SomaticClinicalImpact/ConditionList \
                -block TraitSet -deq '\n' -def 'NA'   -CONTRIB TraitSet@ContributesToAggregateClassification   -EVIDEN TraitSet@LowerLevelOfEvidence -TSYPE TraitSet@Type \
                    -section Trait -deq '\n' -def 'NA'  -element  '&KEYVCV'  '&TSID'  '&TSTYPE' Trait@ID Trait@Type '&CONTRIB' \
                        -subset Trait/XRef -if XRef@DB -equals 'MedGen'  -element  XRef@ID  \
            -group OncogenicityClassification/ConditionList \
                -block TraitSet -deq '\n' -def 'NA'   -CONTRIB TraitSet@ContributesToAggregateClassification   -EVIDEN TraitSet@LowerLevelOfEvidence -TSYPE TraitSet@Type \
                    -section Trait -deq '\n' -def 'NA'  -element  '&KEYVCV'  '&TSID'  '&TSTYPE' Trait@ID Trait@Type '&CONTRIB' \
                        -subset Trait/XRef -if XRef@DB -equals 'MedGen'  -element  XRef@ID    > ~{GENE_NAME}_traitset.txt


    >>>

    output {
        File traitset  = "~{GENE_NAME}_traitset.txt"
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




task extract_clinvar_variants_basic {
    input {
        String GENE_NAME
        File basicxml  
        Int memSizeGB = 10
        Int threadCount = 1
        Int diskSizeGB = 5*round(size(basicxml, "GB")) + 15

    }

    command <<<
        set -eux -o pipefail

        cat  ~{basicxml} |
        xtract -pattern VariationArchive -def "NA" -KEYVCV VariationArchive@Accession -KEYCHANGE "(unknown)" -KEYCONS "(unknown)" -KEYVNAME VariationArchive@VariationName -KEYVDC VariationArchive@DateCreated -KEYVDLU VariationArchive@DateLastUpdated -KEYVMRS VariationArchive@MostRecentSubmission -KEYVTYPE VariationArchive@VariationType -lbl "VCV" -element VariationArchive@Accession VariationArchive@VariationName VariationArchive@VariationType VariationArchive@NumberOfSubmissions VariationArchive@Version \
            -group ClassifiedRecord/SimpleAllele/Location/SequenceLocation  -if SequenceLocation@forDisplay -equals true -def "NA" \
                -KEYASM SequenceLocation@Assembly -KEYCHR SequenceLocation@Chr  -KEYSTART SequenceLocation@start -KEYSTOP SequenceLocation@stop -KEYVCF SequenceLocation@positionVCF -KEYREFA SequenceLocation@referenceAlleleVCF -KEYALTA SequenceLocation@alternateAlleleVCF -KEYVLEN SequenceLocation@variantLength \
                -element SequenceLocation@Assembly SequenceLocation@Chr SequenceLocation@start SequenceLocation@stop SequenceLocation@positionVCF SequenceLocation@referenceAlleleVCF SequenceLocation@alternateAlleleVCF SequenceLocation@variantLength \
            -group ClassifiedRecord/SimpleAllele -element "&KEYVDC"  "&KEYVDLU"  "&KEYVMRS"  \
            -group ClassifiedRecord/Classifications -if GermlineClassification -def "NA" -element GermlineClassification/ReviewStatus GermlineClassification/Description  \
                    -else -lbl "NA\tNA" \
            -group ClassifiedRecord/Classifications -if OncogenicityClassification -def "NA" -element OncogenicityClassification/ReviewStatus OncogenicityClassification/Description  \
                    -else -lbl "NA\tNA" \
            -group ClassifiedRecord/Classifications -if SomaticClinicalImpact -def "NA" -element SomaticClinicalImpact/ReviewStatus SomaticClinicalImpact/Description  \
                    -else -lbl "NA\tNA" \
            -group ClassifiedRecord/SimpleAllele/HGVSlist/HGVS -if NucleotideExpression@MANESelect -equals true \
                -def "NA" -KEYCHANGE NucleotideExpression@change  -KEYCONS -first MolecularConsequence@Type \
            -group ClassifiedRecord/SimpleAllele -element "&KEYCONS"  "&KEYCHANGE"  \
            -group ClassifiedRecord/ClinicalAssertionList/ClinicalAssertion   \
                -deq "\n" -def "NA" -lbl "SCV" -element "&KEYVCV" "&KEYVNAME" "&KEYVTYPE" ClinicalAssertion/ClinVarAccession@Accession ClinicalAssertion/ClinVarAccession@Version "&KEYASM" "&KEYCHR" "&KEYSTART" "&KEYSTOP" "&KEYVCF" "&KEYREFA" "&KEYALTA" "&KEYVLEN" \
                ClinicalAssertion@DateCreated ClinicalAssertion@DateLastUpdated ClinicalAssertion@SubmissionDate \
                Classification/ReviewStatus Classification/GermlineClassification  \
                Classification/ReviewStatus Classification/OncogenicityClassification \
                Classification/ReviewStatus Classification/SomaticClinicalImpact \
                "&KEYCONS"   "&KEYCHANGE" \
                Classification/Comment FunctionalConsequence@Value FunctionalConsequence/Comment  ClinicalAssertion@ID \
                    -block ObservedInList -def "NA"  \
                        -subset ObservedIn/Method -if ObsMethodAttribute/Attribute@Type -equals MethodResult -first ObsMethodAttribute/Attribute |
        sed "s/&gt;/>/g" |
        sed "s/&lt;/</g" | 
        sed "s/â€™/'/g" | 
        sed "s/â€˜/'/g" | 
        sed "s/&amp;/&/g" > ~{GENE_NAME}_basic_res.txt


    >>>

    output {
        File basiccv  = "~{GENE_NAME}_basic_res.txt"
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




task merge_clinvar_variants {
    input {
        Int memSizeGB = 6
        Int threadCount = 1
        File basiccv  
        File traitmap
        File traitset
        Int diskSizeGB = 5*round(size(basiccv, "GB") + size(traitmap, 'GB') + size(traitset, 'GB')) + 2

    }

    command <<<
        set -eux -o pipefail

        python3 <<CODE
        import pandas as pd

        Gene_CV_basic = pd.read_csv("~{basiccv}", delimiter="\t", 
                                    names =["row_type", "VCV_ID", "ClinVar_variant_ID", "variant_type","submissions", "version", 
                                        "assembly", "Chr", "start", "stop",  "pos_VCF", "ref", "alt","variant_length",  
                                        "date_created", "date_updated",  "date_submitted", 
                                        "review_status","germline_classification",
                                        "onco_review_status","oncogenicity_classification",
                                        "som_review_status","somatic_classification",
                                        "variant_effect","txpt_hgvsc",
                                        "comment",
                                        "functional_category", "FA_comment",
                                        "CA_ID", "functional_result"] , header=None, keep_default_na=False)

        trait_set = pd.read_csv("~{traitset}", delimiter="\t", 
                            names = ["VCV_ID", "TraitSet_ID","TS_Type" , "Trait_ID","Trait_Type",  "ContributesToAggregateClassification","MG_ID" ] )
        trait_set['TraitSet_ID'] = trait_set['TraitSet_ID'].fillna("none")
        trait_set[['TraitSet_ID' ]] = trait_set[['TraitSet_ID' ]].astype('str')
        trait_set[['Trait_ID' ]] = trait_set[['Trait_ID' ]].astype('str')
        trait_set['MG_ID'] = trait_set['MG_ID'].fillna("none")
        trait_set['Trait_Type'] = trait_set['Trait_Type'].fillna("none")
        trait_set['TS_Type'] = trait_set['TS_Type'].fillna("none")
        trait_set=trait_set.sort_values([ "VCV_ID", "MG_ID"])
        trait_set['MG_ID'] = trait_set[["VCV_ID", "TraitSet_ID","TS_Type" , "Trait_ID","Trait_Type",  "ContributesToAggregateClassification","MG_ID"]].groupby(["VCV_ID", "TraitSet_ID" ])['MG_ID'].transform(lambda x: '_'.join(x))
        trait_set['Trait_ID'] = trait_set[["VCV_ID", "TraitSet_ID","TS_Type" , "Trait_ID","Trait_Type",  "ContributesToAggregateClassification","MG_ID"]].groupby(["VCV_ID", "TraitSet_ID" ])['Trait_ID'].transform(lambda x: '&'.join(x))
        trait_set =  trait_set[["VCV_ID", "TraitSet_ID","TS_Type" , "Trait_ID","Trait_Type",  "ContributesToAggregateClassification","MG_ID"]].drop_duplicates()


        trait_map = pd.read_csv("~{traitmap}", delimiter="\t", 
                            names = ["label", "VCV_ID", "CA_ID", "Trait_Type", "MG_ID", "MG_disease_name"] )
        trait_map['MG_ID'] = trait_map['MG_ID'].fillna("none")
        trait_map['Trait_Type'] = trait_map['Trait_Type'].fillna("none")
        trait_map['MG_disease_name'] = trait_map['MG_disease_name'].fillna("none")
        trait_map[['CA_ID' ]] = trait_map[['CA_ID' ]].astype('str')
        trait_map = trait_map.sort_values([ "VCV_ID", "MG_ID"])
        trait_map['MG_ID'] = trait_map[["label", "VCV_ID", "CA_ID", "Trait_Type", "MG_ID", "MG_disease_name"]].groupby(["label", "VCV_ID","Trait_Type", "CA_ID"])['MG_ID'].transform(lambda x: '_'.join(x))
        trait_map['MG_disease_name'] = trait_map[["label", "VCV_ID", "CA_ID", "Trait_Type", "MG_ID", "MG_disease_name"]].groupby(["label", "VCV_ID","Trait_Type", "CA_ID","MG_ID"])['MG_disease_name'].transform(lambda x: '&'.join(x))
        trait_map = trait_map[["label", "VCV_ID", "CA_ID", "Trait_Type", "MG_ID", "MG_disease_name"]].drop_duplicates()


        trait_comb = pd.merge(trait_set, trait_map, how='outer', on=["VCV_ID", "MG_ID"])
        trait_comb['MG_ID'] = trait_comb['MG_ID'].fillna("none")
        trait_comb['CA_ID'] = trait_comb['CA_ID'].fillna("none")
        trait_comb['MG_disease_name'] = trait_comb['MG_disease_name'].fillna("none")
        trait_comb['TraitSet_ID'] = trait_comb['TraitSet_ID'].fillna("none")
        trait_comb['Trait_Type_y'] = trait_comb['Trait_Type_y'].fillna("none")
        trait_comb['Trait_Type_x'] = trait_comb['Trait_Type_x'].fillna("none")
        trait_comb['TS_Type'] = trait_comb['TS_Type'].fillna("none")
        trait_comb = trait_comb[trait_comb.CA_ID != "none"]
        trait_comb = trait_comb.drop(trait_comb[(trait_comb['TraitSet_ID'] == "none") & (trait_comb['Trait_Type_y'] == "Finding")].index)

        clinvar_complete = pd.merge(Gene_CV_basic, trait_comb, how='outer', on=["VCV_ID", "CA_ID"])
        clinvar_complete['Chr'] = 'chr' + clinvar_complete['Chr'].astype(str)
        clinvar_complete['CERFAC_variant_id_orig'] = clinvar_complete[['assembly', 'Chr','start','stop','ref','alt' ]].astype(str).agg(':'.join, axis=1)
        clinvar_complete['CERFAC_variant_id_VCF'] = clinvar_complete[['assembly', 'Chr','pos_VCF','ref','alt' ]].astype(str).agg(':'.join, axis=1)
        clinvar_complete['CERFAC_variant_id_HGVS_long'] = clinvar_complete[['assembly', 'Chr','ClinVar_variant_ID' ]].astype(str).agg(':'.join, axis=1)
        clinvar_complete['CERFAC_variant_id_HGVS_short'] = clinvar_complete[['assembly', 'Chr','pos_VCF','txpt_hgvsc' ]].astype(str).agg(':'.join, axis=1)

        cols = ['row_type',
        'VCV_ID',
        'CERFAC_variant_id_orig',
        'CERFAC_variant_id_VCF',
        'CERFAC_variant_id_HGVS_long',
        'CERFAC_variant_id_HGVS_short',
        'ClinVar_variant_ID','submissions',
        #'assembly','Chr',
        'start','stop','pos_VCF','ref','alt',
        'variant_length', 'MG_disease_name','ContributesToAggregateClassification',
        'variant_type', 'variant_effect','txpt_hgvsc', 'review_status',
        'germline_classification', 'oncogenicity_classification',
        'somatic_classification',
        'comment',
        'functional_category',
        'FA_comment',
        'functional_result',
        'date_created',
        'date_updated',
        'date_submitted',
        'version']

        clinvar_complete = clinvar_complete[cols]
        clinvar_complete.to_csv("clinvar_variants.csv", sep=',', index=False )




        CODE

    >>>

    output {
        File clinvar_var = "clinvar_variants.csv"
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






task merge_variants {
    input {
        File gnomadvar
        File clinvar_var
        Int memSizeGB = 6
        Int threadCount = 1
        Int diskSizeGB = 5*round(size(gnomadvar, "GB")) + 20
        String GENE_NAME
    }

    command <<<
        set -eux -o pipefail
        python3 <<CODE
        # python3 ~{MERGE_ALL_SCRIPT} -g ~{gnomadvar} -c ~{clinvar_var}  -o ~{GENE_NAME}_combined_variants.csv
        import os
        import io
        import json
        import string
        import pandas as pd



        cv_table = pd.read_csv("~{clinvar_var}", sep=',' )
        gnomad_vars = pd.read_csv("~{gnomadvar}", sep=',' )
        combined = gnomad_vars.set_index('txpt_hgvsc').join(cv_table.set_index('txpt_hgvsc'), how='outer', lsuffix='_gnomad', rsuffix='_clinvar' )
        combined.sort_values(['txpt_hgvsc'])
        combined.to_csv( "~{GENE_NAME}_combined_variants.csv",  sep=',', index=True )



        CODE

    >>>

    output {
        File combined_var = "~{GENE_NAME}_combined_variants.csv"
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




