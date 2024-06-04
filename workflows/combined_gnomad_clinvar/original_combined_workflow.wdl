version 1.0

workflow annotate_functional_variants {

    meta {
        author: "Allison Cheney"
        email: "archeney@ucsc.edu"
        description: "Extract missense Clinvar and gnomad variants for a specified gene"
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
        Boolean? ETHNICITY_2_GROUP = true
    }
    #The order in which the workflow block and task definitions are arranged in the script does not matter. 
    #Nor does the order of the call statements matter, as we'll see further on.
    call get_gene_info {
        input: GENE_NAME=GENE_NAME

    }
    call get_clinvar_variants{
        input: GENE_NAME=GENE_NAME
    }
    call merge_clinvar_variants{
        input: GENE_NAME=GENE_NAME
    }
    call get_gnomad_variants{
        input: 
            GENE_NAME=GENE_NAME, 
            CHR_ID=get_gene_info.CHR_ID, 
            GENE_START_LOCUS=get_gene_info.GENE_START_LOCUS, 
            GENE_END_LOCUS=get_gene_info.GENE_END_LOCUS

    }
    call merge_variants{
        input:
            gnomadvar=get_gnomad_variants.gnomadvar,
            clinvar_var=merge_clinvar_variants.clinvar_var
    }


    output{
        File combined_var = merge_variants.combined_var
    }

}




task get_gene_info {
    input {
        String GENE_NAME
    }

    command <<<
        esearch -db gene -query "~{GENE_NAME} [GENE] AND homo sapiens [ORGN]" | efetch -format native -mode xml |
        xtract -pattern Entrezgene-Set  -def "-" \
            -element SubSource/SubSource_subtype@value SubSource/SubSource_name \
                -group Entrezgene_locus/Gene-commentary \
                    -if Gene-commentary_heading  -contains "GRCh38"  -def "-" -element Gene-commentary_heading Seq-interval_from Seq-interval_to > gene_info.txt
        echo "getting chromsome of $GENE_NAME"
        # some way to output number of assemblies
        awk '{ print NR-1}' gene_info.txt
        #the chromosome number will always be the first line...
        #echo $CHR_ID
        #CHR_ID=$( awk -F '\t' 'NR == 1 {print $2}' gene_info.txt)
        awk -F '\t' 'NR == 1 {print $2}' gene_info.txt | tee CHR_ID
        
        #... but the remaining variables may have multiple assemblies,
        # and will be printed on multiple lines in an unknown order
        # so we will have to be more specific with NR
        if grep -q -m 1 "GRCh38" gene_info.txt; then
            ASSEMBLY='GRCh38' #
            grep 'GRCh38' gene_info.txt | awk -F '\t' '{print $4}' | tee GENE_START_LOCUS
            grep 'GRCh38' gene_info.txt | awk -F '\t' '{print $5}' | tee GENE_END_LOCUS
        else 
            echo "couldn't find hg38 reference"
        fi
        pwd
    >>>

    output {
        File gene_results = "gene_info.txt"
        String CHR_ID = read_string("CHR_ID")
        Int GENE_START_LOCUS = read_int("GENE_START_LOCUS")
        Int GENE_END_LOCUS = read_int("GENE_END_LOCUS")
        #String ASSEMBLY = 'GRCh38'
        }
    runtime {
        memory: memSizeGB + " GB"
        docker: "ubuntu:18.04"
        preemptible: 1
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

        esearch -db clinvar -query "~{GENE_NAME}[gene] AND single_gene [PROP] AND homo sapiens [ORGN] AND var single nucleotide [FILT]" |
        efetch -format variationid | 
        xtract -pattern VariationArchive -def "NA" -KEYVCV VariationArchive@Accession  -KEYVNAME VariationArchive@VariationName -KEYVDC VariationArchive@DateCreated -KEYVDLU VariationArchive@DateLastUpdated -KEYVMRS VariationArchive@MostRecentSubmission -KEYVTYPE VariationArchive@VariationType -lbl "VCV" -element VariationArchive@Accession VariationArchive@VariationName VariationArchive@VariationType VariationArchive@NumberOfSubmissions VariationArchive@Version \
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
        efetch -format variationid | 
        xtract -pattern RCVAccession -element RCVAccession@Accession > ~{GENE_NAME}_rcv_list.txt

        efetch -db clinvar -input clinvar_rcv_list.txt -format clinvarset | 
        xtract -pattern ClinVarSet -def "NA" -KEYVCV MeasureSet@Acc -KEYRCV  ReferenceClinVarAssertion/ClinVarAccession@Acc \
            -block ClinVarAssertion -deq "\n" -def "NA" -element  "&KEYVCV"  "&KEYRCV" ClinVarAssertion@ID ClinVarAccession@Acc  > ~{GENE_NAME}_rcv_CAID.txt 

        efetch -db clinvar -input clinvar_rcv_list.txt -format clinvarset | 
        xtract -pattern ClinVarSet -def "NA" -KEYVCV MeasureSet@Acc -KEYRCV  ReferenceClinVarAssertion/ClinVarAccession@Acc \
            -group ReferenceClinVarAssertion/TraitSet   -TSID TraitSet@ID \
                -block TraitSet/Trait -deq "\n" -def "NA"  -element  "&KEYVCV"  "&KEYRCV" "&TSID" Trait@ID  \
                    -subset Trait/Name -if ElementValue@Type -equals "Preferred"   -def "NA"  -element  ElementValue \
                    -subset Trait/XRef -if XRef@DB -equals "MedGen"   -def "NA"  -element  XRef@ID  > ~{GENE_NAME}_rcv_trait.txt


    >>>

    output {
        File basiccv  = "~{GENE_NAME}_basic_res.txt"
        File traitset  = "~{GENE_NAME}_traitset.txt"
        File rcvlist  = "~{GENE_NAME}_rcv_list.txt"
        File rcv_CAID  = "~{GENE_NAME}_rcv_CAID.txt"
        File rcv_trait  = "~{GENE_NAME}_rcv_trait.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "ubuntu:18.04"
        preemptible: 1
    }
}





task merge_clinvar_variants {
    input {
        Int memSizeGB = 4
        Int threadCount = 1
        Int diskSizeGB = 5*round(size(input_vcf, "GB")) + 20
        File basiccv  
        File traitmapping 
        File traitset  
        String GENE_NAME

    }

    String basen = sub(sub(basename(input_vcf), ".vcf.bgz$", ""), ".vcf.gz$", "")

    command <<<
        set -eux -o pipefail

        python3 cv_merge_script.py -f ~{basiccv} -t ~{traitmapping} -s ~{traitset}  -o  clinvar_variants.csv


    >>>

    output {
        File clinvar_var = "clinvar_variants.csv"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "ubuntu:18.04"
        preemptible: 1
    }
}




task get_gnomad_variants {
    input {
        Int memSizeGB = 4
        Int threadCount = 1
        Int diskSizeGB = 5*round(size(input_vcf, "GB")) + 20
        String GENE_NAME
        String CHR_ID
        Int GENE_START_LOCUS
        Int GENE_END_LOCUS
    }

    command <<<
        set -eux -o pipefail

        python3 gnomad_variants_script.py -g ~{GENE_NAME} -c ~{CHR_ID} -b ~{GENE_START_LOCUS} -e ~{GENE_END_LOCUS} -o gnomad_variants_MANE.csv
    >>>

    output {
        File gnomadvar = "gnomad_variants_MANE.csv"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "ubuntu:18.04"
        preemptible: 1
    }
}



task merge_variants {
    input {
        File gnomadvar
        File clinvar_var
        Int memSizeGB = 4
        Int threadCount = 1
        Int diskSizeGB = 5*round(size(input_vcf, "GB")) + 20
        String GENE_NAME
        String CHR_ID
        Int GENE_START_LOCUS
        Int GENE_END_LOCUS
    }

    command <<<
        set -eux -o pipefail

        python3 merge_all_variants.py -g ~{GENE_NAME} -o ~{GENE_NAME}_combined_variants.csv
    >>>

    output {
        File combined_var = "~{GENE_NAME}_combined_variants.csv"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "ubuntu:18.04"
        preemptible: 1
    }
}




