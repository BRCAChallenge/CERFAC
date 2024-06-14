version 1.0

workflow annotate_functional_variants {

    meta {
        author: "Allison Cheney"
        email: "archeney@ucsc.edu"
        description: "Extract missense Clinvar  variants for a specified gene"
    }

    parameter_meta {
        GENE_NAME: "Name of relevant gene"

    }
    input {
        String GENE_NAME    
    }
    #The order in which the workflow block and task definitions are arranged in the script does not matter. 
    #Nor does the order of the call statements matter, as we'll see further on.
    call get_clinvar_variants_file{
        input: GENE_NAME=GENE_NAME
    }
    call extract_clinvar_variants_basic{
        input: 
            GENE_NAME=GENE_NAME,
            basicxml=get_clinvar_variants_file.basicxml
    }
    call extract_clinvar_variants_traitset{
        input: 
            GENE_NAME=GENE_NAME,
            basicxml=get_clinvar_variants_file.basicxml
    }
    call extract_clinvar_variants_traitmap{
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
        File clinvar_var = merge_clinvar_variants.clinvar_var
    }

}


task get_clinvar_variants_file {
    input {
        String GENE_NAME
        Int memSizeGB = 6
        Int threadCount = 1
        Int diskSizeGB = 25
    }

    command <<<
        set -eux -o pipefail

        esearch -db clinvar -query "~{GENE_NAME}[GENE] AND single_gene [PROP] AND homo sapiens [ORGN] AND var single nucleotide [FILT]" | efetch -format variationid -start 1 -stop 1 | 
        xtract -pattern VariationArchive  \
            -group ClassifiedRecord/SimpleAllele/GeneList/Gene/Location/SequenceLocation  -if SequenceLocation@Assembly -equals "GRCh38" -def "NA" \
                -element SequenceLocation@Assembly  SequenceLocation@Chr SequenceLocation@start  SequenceLocation@stop > ~{GENE_NAME}_positions.txt

        esearch -db clinvar -query "~{GENE_NAME}[GENE] AND single_gene [PROP] AND homo sapiens [ORGN] AND var single nucleotide [FILT]" |
        efetch -format variationid > basic.xml
    >>>

    output {
        File basicxml  = "basic.xml"
        File gene_positions = "~{GENE_NAME}_positions.txt"
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
        Int memSizeGB = 6
        Int threadCount = 1
        Int diskSizeGB = 5*round(size(basicxml, "GB")) + 2

    }

    command <<<
        set -eux -o pipefail

        xtract -input ~{basicxml} -pattern VariationArchive -def "NA" -KEYVCV VariationArchive@Accession  -KEYVNAME VariationArchive@VariationName -KEYVDC VariationArchive@DateCreated -KEYVDLU VariationArchive@DateLastUpdated -KEYVMRS VariationArchive@MostRecentSubmission -KEYVTYPE VariationArchive@VariationType -lbl "VCV" -element VariationArchive@Accession VariationArchive@VariationName VariationArchive@VariationType VariationArchive@NumberOfSubmissions VariationArchive@Version \
             -group ClassifiedRecord/SimpleAllele/GeneList/Gene/Location/SequenceLocation  -if SequenceLocation@Assembly -equals "GRCh38" -def "NA" \
                -element SequenceLocation@Assembly SequenceLocation@Chr SequenceLocation@start SequenceLocation@stop \
            -group ClassifiedRecord/SimpleAllele/Location/SequenceLocation  -if SequenceLocation@forDisplay -equals true -def "NA" \
                -KEYASM SequenceLocation@Assembly -KEYCHR SequenceLocation@Chr  -KEYSTART SequenceLocation@start -KEYSTOP SequenceLocation@stop -KEYREFA SequenceLocation@referenceAlleleVCF -KEYALTA SequenceLocation@alternateAlleleVCF -KEYVLEN SequenceLocation@variantLength \
                -element SequenceLocation@Assembly SequenceLocation@Chr SequenceLocation@start SequenceLocation@stop SequenceLocation@referenceAlleleVCF SequenceLocation@alternateAlleleVCF SequenceLocation@variantLength \
            -group ClassifiedRecord/SimpleAllele/HGVSlist/HGVS -if NucleotideExpression@MANESelect -equals true \
                -def "NA" -KEYCHANGE NucleotideExpression@change  -KEYCONS -first MolecularConsequence@Type -first MolecularConsequence@Type NucleotideExpression@change   \
            -group ClassifiedRecord/SimpleAllele/ -element "&KEYVDC"  "&KEYVDLU"  "&KEYVMRS"  \
            -group ClassifiedRecord/Classifications -if GermlineClassification -def "NA" -element GermlineClassification/ReviewStatus GermlineClassification/Description  \
                    -else -lbl "NA\tNA" \
            -group ClassifiedRecord/Classifications -if OncogenicityClassification -def "NA" -element OncogenicityClassification/ReviewStatus OncogenicityClassification/Description  \
                    -else -lbl "NA\tNA" \
            -group ClassifiedRecord/Classifications -if SomaticClinicalImpact -def "NA" -element SomaticClinicalImpact/ReviewStatus SomaticClinicalImpact/Description  \
                    -else -lbl "NA\tNA" \
            -group ClassifiedRecord/ClinicalAssertionList/ClinicalAssertion   \
                -deq "\n" -def "NA" -lbl "SCV" -element "&KEYVCV" "&KEYVNAME" "&KEYVTYPE" ClinicalAssertion/ClinVarAccession@Accession ClinicalAssertion/ClinVarAccession@Version "&KEYASM" "&KEYCHR" "&KEYSTART" "&KEYSTOP"  "&KEYREFA" "&KEYALTA" "&KEYVLEN" "&KEYCONS"   "&KEYCHANGE"\
                ClinicalAssertion@DateCreated ClinicalAssertion@DateLastUpdated ClinicalAssertion@SubmissionDate \
                Classification/ReviewStatus Classification/GermlineClassification  \
                Classification/ReviewStatus Classification/OncogenicityClassification \
                Classification/ReviewStatus Classification/SomaticClinicalImpact @ClinicalImpactAssertionType @ClinicalImpactClinicalSignificance @DrugForTherapeuticAssertion \
                Classification/Comment FunctionalConsequence@Value FunctionalConsequence/Comment ClinVarAccession@OrganizationCategory Classification/StudyDescription ClinicalAssertion@ID \
                    -block ObservedInList -def "NA" -first Method/Description Method/MethodType Sample/CellLine \
                        -subset ObservedIn/Method -if ObsMethodAttribute/Attribute@Type -equals MethodResult -first ObsMethodAttribute/Attribute > basicintmd.txt
        
        sed "s/&gt;/>/g" basicintmd.txt | 
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
        maxRetries: 4
    }
}


task extract_clinvar_variants_traitset {
    input {
        String GENE_NAME
        File basicxml  
        Int memSizeGB = 6
        Int threadCount = 1
        Int diskSizeGB = 5*round(size(basicxml, "GB")) + 2

    }

    command <<<
        set -eux -o pipefail

        xtract -input ~{basicxml} -pattern VariationArchive -def 'NA' -KEYVCV VariationArchive@Accession \
            -group GermlineClassification/ConditionList \
                -block TraitSet -deq '\n' -def 'NA'   -TSID TraitSet@ID  -CONTRIB TraitSet@ContributesToAggregateClassification    \
                    -subset Trait -deq '\n' -element  '&KEYVCV'  '&TSID' Trait@ID '&CONTRIB' \
            -group SomaticClinicalImpact/ConditionList \
                -block TraitSet -deq '\n' -def 'NA'  -TSID TraitSet@ID  -CONTRIB TraitSet@ContributesToAggregateClassification   -EVIDEN TraitSet@LowerLevelOfEvidence -TTYPE TraitSet@Type \
                    -subset Trait -deq '\n' -def 'NA' -element  '&KEYVCV'  '&TSID' Trait@ID '&CONTRIB' '&EVIDEN' '&TTYPE'\
            -group OncogenicityClassification/ConditionList \
                -block TraitSet -deq '\n' -def 'NA'   -TSID TraitSet@ID  -CONTRIB TraitSet@ContributesToAggregateClassification   -EVIDEN TraitSet@LowerLevelOfEvidence -TTYPE TraitSet@Type \
                    -subset Trait -deq '\n' -def 'NA'  -element  '&KEYVCV'  '&TSID' Trait@ID '&CONTRIB' '&EVIDEN' '&TTYPE'  > ~{GENE_NAME}_traitset.txt


    >>>

    output {
        File traitset  = "~{GENE_NAME}_traitset.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "allisoncheney/cerfac_terra:clinvar"
        maxRetries: 4
    }
}

task extract_clinvar_variants_traitmap {
    input {
        String GENE_NAME
        File basicxml  
        Int memSizeGB = 6
        Int threadCount = 1
        Int diskSizeGB = 5*round(size(basicxml, "GB")) + 2

    }

    command <<<
        set -eux -o pipefail

        xtract -input~{basicxml} -pattern VariationArchive -def 'NA' -KEYVCV VariationArchive@Accession \
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
        maxRetries: 4
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

        python3 cv_merge_script.py -f ~{basiccv} -m ~{traitmap} -s ~{traitset}  -o  clinvar_variants.csv


    >>>

    output {
        File clinvar_var = "clinvar_variants.csv"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "allisoncheney/cerfac_terra:clinvar"
        maxRetries: 4
    }
}
