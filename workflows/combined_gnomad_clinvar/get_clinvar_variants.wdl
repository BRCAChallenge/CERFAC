version 1.0

workflow get_clinvar_variants {

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

    output{
        File output_clinvar_variants = merge_clinvar_variants.clinvar_var
        String output_clinvar_variants_count = merge_clinvar_variants.clinvar_variants_count
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

        esearch -db clinvar -query "~{GENE_NAME}[GENE] AND homo sapiens [ORGN] AND (varlen 49 or less[FILTER]) NOT (near gene upstream[PROP]) NOT (near gene downstream[PROP])" |
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
        maxRetries: 1
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
                -group TraitMapping -deq '\n' -def 'None given' -lbl 'traitmapping' -element '&KEYVCV' @ClinicalAssertionID @TraitType MedGen@CUI MedGen@Name > ~{GENE_NAME}_traitmapping.txt

    >>>

    output {
        File traitmap  = "~{GENE_NAME}_traitmapping.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "allisoncheney/cerfac_terra:clinvar"
        maxRetries: 1
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
        xtract -pattern VariationArchive -def "NA" -KEYVCV VariationArchive@Accession -KEYCHANGE "(unknown)" -KEYCONS "(unknown)" -KEYVNAME VariationArchive@VariationName -KEYVDC VariationArchive@DateCreated -KEYVDLU VariationArchive@DateLastUpdated -KEYVTYPE VariationArchive@VariationType -KEYSUBNUM VariationArchive@NumberOfSubmissions \
        -KEYGREVOV "(None given)" -KEYGCLASSOV "(None given)" -KEYOREVOV "(None given)" -KEYOCLASSOV "(None given)" -KEYSREVOV "(None given)" -KEYSCLASSOV "(None given)"\
            -group ClassifiedRecord/SimpleAllele/Location/SequenceLocation  -if SequenceLocation@forDisplay -equals true -def "NA" \
                -KEYASM SequenceLocation@Assembly -KEYCHR SequenceLocation@Chr  -KEYSTART SequenceLocation@start -KEYSTOP SequenceLocation@stop -KEYVCF SequenceLocation@positionVCF -KEYREFA SequenceLocation@referenceAlleleVCF -KEYALTA SequenceLocation@alternateAlleleVCF -KEYVLEN SequenceLocation@variantLength \
            -group ClassifiedRecord/Classifications -if GermlineClassification  -KEYGREVOV GermlineClassification/ReviewStatus -KEYGCLASSOV GermlineClassification/Description  \
            -group ClassifiedRecord/Classifications -if OncogenicityClassification  -KEYOREVOV OncogenicityClassification/ReviewStatus -KEYOCLASSOV OncogenicityClassification/Description  \
            -group ClassifiedRecord/Classifications -if SomaticClinicalImpact  -KEYSREVOV SomaticClinicalImpact/ReviewStatus -KEYSCLASSOV SomaticClinicalImpact/Description  \
            -group ClassifiedRecord/SimpleAllele/HGVSlist/HGVS -if NucleotideExpression@MANESelect -equals true \
                 -KEYCHANGE NucleotideExpression@change  -KEYCONS -first MolecularConsequence@Type \
            -group ClassifiedRecord/ClinicalAssertionList/ClinicalAssertion   \
                -deq "\n" -def "None given" -element "&KEYVCV" "&KEYVNAME" "&KEYVTYPE" "&KEYSUBNUM" ClinicalAssertion/ClinVarAccession@Accession "&KEYASM" "&KEYCHR" "&KEYSTART" "&KEYSTOP" "&KEYVCF" "&KEYREFA" "&KEYALTA" "&KEYVLEN" \
                "&KEYVDC"  "&KEYVDLU" ClinicalAssertion@DateCreated ClinicalAssertion@DateLastUpdated ClinicalAssertion@SubmissionDate \
                "&KEYGREVOV" "&KEYGCLASSOV" "&KEYOREVOV" "&KEYOCLASSOV" "&KEYSREVOV" "&KEYSCLASSOV" Classification/ReviewStatus Classification/GermlineClassification  \
                Classification/OncogenicityClassification Classification/SomaticClinicalImpact \
                "&KEYCONS"   "&KEYCHANGE" \
                Classification/Comment FunctionalConsequence@Value FunctionalConsequence/Comment  ClinicalAssertion@ID \
                    -block ObservedInList -def "None given"  \
                        -subset ObservedIn/Method -if ObsMethodAttribute/Attribute@Type -equals MethodResult -first ObsMethodAttribute/Attribute  |
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
        maxRetries: 1
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
        Int diskSizeGB = 5*round(size(basiccv, "GB") + size(traitmap, 'GB') + size(traitset, 'GB')) + 5

    }
    #if you delete any columns make sure to change the number of columns included
    command <<<
        set -eux -o pipefail

        python3 <<CODE
        import pandas as pd

        Gene_CV_basic = pd.read_csv("~{basiccv}", delimiter="\t", engine='python',
                                    names =["VCV_ID", "ClinVar_variant_ID", "variant_class","number_submissions", "SCV_ID", 
                                        "assembly", "Chr", "start", "stop",  "pos_VCF", "ref", "alt","variant_length",  
                                        "date_variant_created", "date_variant_updated",  "date_submission_created", "date_submission_updated","date_submitted", 
                                        "overall_germline_review_status","overall_germline_classification","overall_onco_review_status","overall_oncogenicity_classification","overall_som_review_status","overall_somatic_classification",
                                        "submission_review_status","submission_germline_classification", "submission_oncogenicity_classification", "submission_somatic_classification",
                                        "variant_effect","txpt_hgvsc",
                                        "comment",
                                        "functional_category", "functional_comment",
                                        "CA_ID", "functional_result", "extra1", "extra2"] ,   header=None, keep_default_na=False)
        Gene_CV_basic.iloc[:, 0:35]

                                

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

        Gene_CV_basic['CA_ID'] = Gene_CV_basic['CA_ID'].fillna("none")
        Gene_CV_basic[['CA_ID' ]] = Gene_CV_basic[['CA_ID' ]].astype('str')
        trait_comb[['CA_ID' ]] = trait_comb[['CA_ID' ]].astype('str')

        clinvar_complete = pd.merge(Gene_CV_basic, trait_comb, how='outer', on=["VCV_ID", "CA_ID"])
        clinvar_complete['Chr'] = 'chr' + clinvar_complete['Chr'].astype(str)
        clinvar_complete['CERFAC_variant_id_VCF'] = clinvar_complete[['assembly', 'Chr','pos_VCF','ref','alt' ]].astype(str).agg(':'.join, axis=1)
        clinvar_complete['CERFAC_variant_id_HGVS_long'] = clinvar_complete[['assembly', 'Chr','ClinVar_variant_ID' ]].astype(str).agg(':'.join, axis=1)
        clinvar_complete['CERFAC_variant_id_HGVS_short'] = clinvar_complete[['assembly', 'Chr','pos_VCF','txpt_hgvsc' ]].astype(str).agg(':'.join, axis=1)
        clinvar_complete['txpt_hgvsc_from_ID'] = clinvar_complete['ClinVar_variant_ID'].str.split(pat=":", n=1,  regex=False).str.get(1)
        clinvar_complete['txpt_hgvsc_from_ID_no_pro'] = clinvar_complete['txpt_hgvsc_from_ID'].str.split(pat=" ", n=1,  regex=False).str.get(0)
        clinvar_complete['hgvs_pro'] = clinvar_complete['ClinVar_variant_ID'].str.split(pat=" ", n=1,  regex=False).str.get(1)
        clinvar_complete['hgvs_pro'] = clinvar_complete['hgvs_pro'].replace(regex=True, to_replace='\)', value='')
        clinvar_complete['hgvs_pro'] = clinvar_complete['hgvs_pro'].replace(regex=True, to_replace='\(', value='')


        cols = ['VCV_ID','txpt_hgvsc_from_ID','hgvs_pro',
        'CERFAC_variant_id_VCF',
        'CERFAC_variant_id_HGVS_long',
        'CERFAC_variant_id_HGVS_short',
        'ClinVar_variant_ID','number_submissions', 'SCV_ID', 
        'start','stop','pos_VCF','ref','alt',
        'variant_length', 'MG_disease_name','ContributesToAggregateClassification',
        'variant_class', 'variant_effect','txpt_hgvsc', 
        'overall_germline_review_status','overall_germline_classification','submission_review_status','submission_germline_classification',
        'overall_onco_review_status','overall_oncogenicity_classification','submission_oncogenicity_classification',
        'overall_som_review_status','overall_somatic_classification','submission_somatic_classification',
        'comment',
        'functional_category','functional_comment','functional_result', 
        'date_variant_created', 'date_variant_updated',  'date_submission_created', 'date_submission_updated', 'date_submitted', 'txpt_hgvsc_from_ID_no_pro']

        clinvar_complete = clinvar_complete[cols]


        clinvar_complete = clinvar_complete[clinvar_complete.variant_effect != "5 prime UTR variant"]
        clinvar_complete = clinvar_complete[clinvar_complete.variant_effect != "3 prime UTR variant"]

        clinvar_variants_count_pd = str(clinvar_complete['txpt_hgvsc'].nunique())
        file_name = "clinvarcount.txt"
        with open(file_name, 'w') as x_file:
            x_file.write(clinvar_variants_count_pd)
        clinvar_complete = clinvar_complete.rename(columns={"ref": "allele_ref",  "alt": "allele_alt",   "start": "pos_start",   "stop": "pos_stop"}, errors='raise')
        clinvar_complete['variant_source']="ClinVar"
        clinvar_complete = clinvar_complete.add_suffix('_clinvar')

        clinvar_complete['txpt_ref_from_ID'] = clinvar_complete['ClinVar_variant_ID_clinvar'].str.split(pat=":", n=1,  regex=False).str.get(0)

        clinvar_complete['ref_txpt_clinvar'] = clinvar_complete['txpt_ref_from_ID'].str.split(pat="(", n=1,  regex=False).str.get(0)
        clinvar_complete['hgvs_cdna_clinvar'] = clinvar_complete[['ref_txpt_clinvar', 'txpt_hgvsc' ]].astype(str).agg(':'.join, axis=1)
        clinvar_complete['hgvs_cdna_clinvar_from_ID'] = clinvar_complete[['ref_txpt_clinvar', 'txpt_hgvsc_from_ID_no_pro_clinvar' ]].astype(str).agg(':'.join, axis=1)


        clinvar_complete.to_csv("clinvar_variants.csv", sep=',', index=False )


        CODE

    >>>

    output {
        File clinvar_var = "clinvar_variants.csv"
        String clinvar_variants_count = read_string("clinvarcount.txt")
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




