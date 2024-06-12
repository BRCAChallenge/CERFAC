version 1.0

workflow annotate_functional_variants {

    meta {
        author: "Allison Cheney"
        email: "archeney@ucsc.edu"
        description: "Extract missense Clinvar  variants for a specified gene"
    }

    parameter_meta {
        GENE_NAME: "Name of relevant gene"
        ASSEMBLY: "Default: GRCh38"
        CHR_ID: "Name of chromosome in number"
        GENE_START_LOCUS: "start position of gene"
        GENE_END_LOCUS: "end position of gene"
        REMOVE_INTRONIC_VAR: "Should variants in introns be removed? Default: true"
        REMOVE_SPLICING_VAR: "Should splicing variants be removed? Default: true"
        REMOVE_NS_ENDPOINT_VAR: "Should variants past the nonsense endpoint be removed? Default: true"
        NS_ENDPOINT: "Point where nonsense variants are not pathogenic. Default: 50 bp before end of penultimate exon"
        ETHNICITY_2_GROUP: "Group ethnicities into Western European vs others. Ashkenazis and Finns excluded. Default: true"

    }
    input {
        
        String GENE_NAME

        Int NS_ENDPOINT = 50

        Boolean? REMOVE_INTRONIC_VAR = true
        Boolean? REMOVE_SPLICING_VAR = true
        Boolean? REMOVE_NS_ENDPOINT_VAR = true
    }
    #The order in which the workflow block and task definitions are arranged in the script does not matter. 
    #Nor does the order of the call statements matter, as we'll see further on.
    call get_clinvar_variants{
        input: GENE_NAME=GENE_NAME
    }
    call merge_clinvar_variants{
        input: 
            GENE_NAME=GENE_NAME,
            basiccv=get_clinvar_variants.basiccv, 
            traitset=get_clinvar_variants.traitset, 
            traitmap=get_clinvar_variants.traitmap
    }

    output{
        File clinvar_var = merge_clinvar_variants.clinvar_var
    }

}


task get_clinvar_variants {
    input {
        String GENE_NAME
        Int memSizeGB = 4
        Int threadCount = 1
        Int diskSizeGB = 5*round(size(input_vcf, "GB")) + 20
    }

    command <<<
        set -eux -o pipefail

        esearch -db clinvar -query "~{GENE_NAME}[gene] AND single_gene [PROP] AND homo sapiens [ORGN] AND var single nucleotide [FILT]" | efetch -format variationid -start 1 -stop 1 | 
        xtract -pattern VariationArchive  \
            -group ClassifiedRecord/SimpleAllele/GeneList/Gene/Location/SequenceLocation  -if SequenceLocation@Assembly -equals "GRCh38" -def "NA" \
                -element SequenceLocation@Assembly  SequenceLocation@Chr SequenceLocation@start  SequenceLocation@stop > ~{GENE_NAME}_positions.txt

        esearch -db clinvar -query "~{GENE_NAME}[gene] AND single_gene [PROP] AND homo sapiens [ORGN] AND var single nucleotide [FILT]" |
        efetch -format variationid | 
        xtract -pattern VariationArchive -def "NA" -KEYVCV VariationArchive@Accession  -KEYVNAME VariationArchive@VariationName -KEYVDC VariationArchive@DateCreated -KEYVDLU VariationArchive@DateLastUpdated -KEYVMRS VariationArchive@MostRecentSubmission -KEYVTYPE VariationArchive@VariationType -lbl "VCV" -element VariationArchive@Accession VariationArchive@VariationName VariationArchive@VariationType VariationArchive@NumberOfSubmissions VariationArchive@Version \
             -group ClassifiedRecord/SimpleAllele/GeneList/Gene/Location/SequenceLocation  -if SequenceLocation@Assembly -equals "GRCh38" -def "NA" \
                -element SequenceLocation@Assembly SequenceLocation@Chr SequenceLocation@start SequenceLocation@stop \
            -group ClassifiedRecord/SimpleAllele/Location/SequenceLocation  -if SequenceLocation@forDisplay -equals true -def "NA" \
                -KEYASM SequenceLocation@Assembly -KEYCHR SequenceLocation@Chr  -KEYSTART SequenceLocation@start -KEYSTOP SequenceLocation@stop -KEYREFA SequenceLocation@referenceAlleleVCF -KEYALTA SequenceLocation@alternateAlleleVCF -KEYVLEN SequenceLocation@variantLength \
                -element SequenceLocation@Assembly SequenceLocation@Chr SequenceLocation@start SequenceLocation@stop SequenceLocation@referenceAlleleVCF SequenceLocation@alternateAlleleVCF SequenceLocation@variantLength \
            -group ClassifiedRecord/SimpleAllele/HGVSlist/HGVS -if NucleotideExpression@MANESelect -equals true \
                -def "NA"  -first MolecularConsequence@Type NucleotideExpression@change  -KEYCHANGE NucleotideExpression@change \
                -KEYCONS -first MolecularConsequence@Type \
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
                        -subset ObservedIn/Method -if ObsMethodAttribute/Attribute@Type -equals MethodResult -first ObsMethodAttribute/Attribute |
        sed "s/&gt;/>/g" | 
        sed "s/&lt;/</g" | 
        sed "s/â€™/'/g" | 
        sed "s/â€˜/'/g" | 
        sed "s/&amp;/&/g" > ~{GENE_NAME}_basic_res.txt

        esearch -db clinvar -query "~{GENE_NAME}[gene] AND single_gene [PROP] AND homo sapiens [ORGN] AND var single nucleotide [FILT]" |    efetch -format variationid | 
        xtract -pattern VariationArchive -def "NA" -KEYVCV VariationArchive@Accession \
            -group GermlineClassification/ConditionList \
                -block TraitSet -deq "\n" -def "NA"   -TSID TraitSet@ID  -CONTRIB TraitSet@ContributesToAggregateClassification    \
                    -subset Trait -deq "\n" -element  "&KEYVCV"  "&TSID" Trait@ID "&CONTRIB" \
            -group SomaticClinicalImpact/ConditionList \
                -block TraitSet -deq "\n" -def "NA"  -TSID TraitSet@ID  -CONTRIB TraitSet@ContributesToAggregateClassification   -EVIDEN TraitSet@LowerLevelOfEvidence -TTYPE TraitSet@Type \
                    -subset Trait -deq "\n" -def "NA" -element  "&KEYVCV"  "&TSID" Trait@ID "&CONTRIB" "&EVIDEN" "&TTYPE"\
            -group OncogenicityClassification/ConditionList \
                -block TraitSet -deq "\n" -def "NA"   -TSID TraitSet@ID  -CONTRIB TraitSet@ContributesToAggregateClassification   -EVIDEN TraitSet@LowerLevelOfEvidence -TTYPE TraitSet@Type \
                    -subset Trait -deq "\n" -def "NA"  -element  "&KEYVCV"  "&TSID" Trait@ID "&CONTRIB" "&EVIDEN" "&TTYPE" |
        sort -k1 -k2  -k3 > ~{GENE_NAME}_traitset.txt

        esearch -db clinvar -query "~{GENE_NAME}[gene] AND single_gene [PROP] AND homo sapiens [ORGN] AND var single nucleotide [FILT]" |
        efetch -format variationid | xtract -pattern VariationArchive -def "NA" -KEYVCV VariationArchive@Accession \
            -group TraitMapping -deq "\n" -def "None" -lbl "traitmapping" -element "&KEYVCV" @ClinicalAssertionID @TraitType MedGen@CUI MedGen@Name |
        sort -k2  -k3 -k4 > ~{GENE_NAME}_traitmapping.txt



    >>>

    output {
        File basiccv  = "~{GENE_NAME}_basic_res.txt"
        File traitset  = "~{GENE_NAME}_traitset.txt"
        File traitmap  = "~{GENE_NAME}_traitmapping.txt"
        File gene_positions = "~{GENE_NAME}_positions.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "allisoncheney/cerfac_terra@sha256:ea88241f4bf4b4e1c8a06a8dc679c9daa32a4dfe5a24ee577cd48d1153291717"
        preemptible: 1
    }
}





task merge_clinvar_variants {
    input {
        Int memSizeGB = 4
        Int threadCount = 1
        Int diskSizeGB = 5*round(size(basiccv, "GB")) + 20
        File basiccv  
        File traitmap
        File traitset  
        String GENE_NAME

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
        docker: "allisoncheney/cerfac_terra@sha256:ea88241f4bf4b4e1c8a06a8dc679c9daa32a4dfe5a24ee577cd48d1153291717"
        preemptible: 1
    }
}
