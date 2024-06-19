version 1.0

workflow annotate_functional_variants {

    meta {
        author: "Allison Cheney"
        email: "archeney@ucsc.edu"
        description: "Extract variants for a specified gene"
    }

    parameter_meta {
        GENE_NAME: "Name of relevant gene"
        GNOMAD_SCRIPT: "gnomad_variants_script.py file used for extracting variants from gnomad"
        MERGE_ALL_SCRIPT: "merge_all_variants.py"
        CVSCRIPT: "cv_merge_script.py"
    }
    input {
        String GENE_NAME
        File CVSCRIPT
        File GNOMAD_SCRIPT
        File MERGE_ALL_SCRIPT
        
    }
    #The order in which the workflow block and task definitions are arranged in the script does not matter. 
    #Nor does the order of the call statements matter, as we'll see further on.
    call extract_gene_loc {
        input: GENE_NAME=GENE_NAME
    }
    call get_gnomad_variants{
        input: 
            GENE_NAME=GENE_NAME, 
            GNOMAD_SCRIPT=GNOMAD_SCRIPT,
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
            CVSCRIPT=CVSCRIPT,
            basiccv=extract_clinvar_variants_basic.basiccv, 
            traitset=extract_clinvar_variants_traitset.traitset, 
            traitmap=extract_clinvar_variants_traitmap.traitmap
    }
    call merge_variants{
        input:
            MERGE_ALL_SCRIPT=MERGE_ALL_SCRIPT,
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
        File GNOMAD_SCRIPT
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

        python3 ~{GNOMAD_SCRIPT} -g ~{GENE_NAME} -c ~{CHR_ID} -b ~{GENE_START_LOCUS} -e ~{GENE_END_LOCUS} -o gnomad_variants_MANE.csv
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

        esearch -db clinvar -query "~{GENE_NAME}[GENE] AND single_gene [PROP] AND homo sapiens [ORGN] AND var single nucleotide [FILT]" |
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
        xtract -pattern VariationArchive -def "NA" -KEYVCV VariationArchive@Accession  -KEYVNAME VariationArchive@VariationName -KEYVDC VariationArchive@DateCreated -KEYVDLU VariationArchive@DateLastUpdated -KEYVMRS VariationArchive@MostRecentSubmission -KEYVTYPE VariationArchive@VariationType -lbl "VCV" -element VariationArchive@Accession VariationArchive@VariationName VariationArchive@VariationType VariationArchive@NumberOfSubmissions VariationArchive@Version \
            -group ClassifiedRecord/SimpleAllele/Location/SequenceLocation  -if SequenceLocation@forDisplay -equals true -def "NA" \
                -KEYASM SequenceLocation@Assembly -KEYCHR SequenceLocation@Chr  -KEYSTART SequenceLocation@start -KEYSTOP SequenceLocation@stop -KEYREFA SequenceLocation@referenceAlleleVCF -KEYALTA SequenceLocation@alternateAlleleVCF -KEYVLEN SequenceLocation@variantLength \
                -element SequenceLocation@Assembly SequenceLocation@Chr SequenceLocation@start SequenceLocation@stop SequenceLocation@referenceAlleleVCF SequenceLocation@alternateAlleleVCF SequenceLocation@variantLength \
            -group ClassifiedRecord/SimpleAllele -element "&KEYVDC"  "&KEYVDLU"  "&KEYVMRS"  \
            -group ClassifiedRecord/Classifications -if GermlineClassification -def "NA" -element GermlineClassification/ReviewStatus GermlineClassification/Description  \
                    -else -lbl "NA\tNA" \
            -group ClassifiedRecord/Classifications -if OncogenicityClassification -def "NA" -element OncogenicityClassification/ReviewStatus OncogenicityClassification/Description  \
                    -else -lbl "NA\tNA" \
            -group ClassifiedRecord/Classifications -if SomaticClinicalImpact -def "NA" -element SomaticClinicalImpact/ReviewStatus SomaticClinicalImpact/Description  \
                    -else -lbl "NA\tNA" \
            -group ClassifiedRecord/SimpleAllele/HGVSlist/HGVS -if NucleotideExpression@MANESelect -equals true \
                -def "NA" -KEYCHANGE NucleotideExpression@change  -KEYCONS -first MolecularConsequence@Type -first MolecularConsequence@Type NucleotideExpression@change \
            -group ClassifiedRecord/ClinicalAssertionList/ClinicalAssertion   \
                -deq "\n" -def "NA" -lbl "SCV" -element "&KEYVCV" "&KEYVNAME" "&KEYVTYPE" ClinicalAssertion/ClinVarAccession@Accession ClinicalAssertion/ClinVarAccession@Version "&KEYASM" "&KEYCHR" "&KEYSTART" "&KEYSTOP"  "&KEYREFA" "&KEYALTA" "&KEYVLEN" \
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
        File CVSCRIPT
        File basiccv  
        File traitmap
        File traitset
        Int diskSizeGB = 5*round(size(basiccv, "GB") + size(traitmap, 'GB') + size(traitset, 'GB')) + 2

    }

    command <<<
        set -eux -o pipefail

        python3 ~{CVSCRIPT} -f ~{basiccv} -m ~{traitmap} -s ~{traitset}  -o  clinvar_variants.csv


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
        File MERGE_ALL_SCRIPT
        File clinvar_var
        Int memSizeGB = 6
        Int threadCount = 1
        Int diskSizeGB = 5*round(size(gnomadvar, "GB")) + 20
        String GENE_NAME
    }

    command <<<
        set -eux -o pipefail

        python3 ~{MERGE_ALL_SCRIPT} -g ~{gnomadvar} -c ~{clinvar_var}  -o ~{GENE_NAME}_combined_variants.csv
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




