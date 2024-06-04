version 1.0

workflow annotate_variants {

    meta {
        author: "Allison Cheney"
        email: "archeney@ucsc.edu"
        description: "Annotate a gene with frequencies in gnomAD and ClinVar"
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
        String ASSEMBLY
        String CHR_ID
        Int GENE_START_LOCUS
        Int GENE_END_LOCUS
        Int NS_ENDPOINT = 50
        Boolean REMOVE_INTRONIC_VAR = true
        Boolean REMOVE_SPLICING_VAR = true
        Boolean REMOVE_NS_ENDPOINT_VAR = true
        Boolean ETHNICITY_2_GROUP = true
    }
    #The order in which the workflow block and task definitions are arranged in the script does not matter. 
    #Nor does the order of the call statements matter, as we'll see further on.

    output{
        File outputfile = exampletask.outputfile
    }

}
task hello_world {
    input {
        File input_vcf
        String GENE_NAME
        Int memSizeGB = 4
    }

  command <<<
    echo "hello world" > result_wdl_test.txt

  >>>

  output {
    File results = "~{sample_id}.out"
  }
    runtime {
        memory: memSizeGB + " GB"
        docker: "ubuntu:18.04"
        preemptible: 1
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
        CHR_ID=$( awk -F '\t' 'NR == 1 {print $2}' gene_info.txt)
        echo $CHR_ID
        #... but the remaining variables may have multiple assemblies,
        # and will be printed on multiple lines in an unknown order
        # so we will have to be more specific with NR
        if grep -q -m 1 "GRCh38" gene_info.txt; then
            ASSEMBLY=GRCh38 #
        else 
            echo "couldn't find hg38 reference"
        fi

        GENE_START_LOCUS=$( grep 'GRCh38' gene_info.txt | awk -F '\t' '{print $4}')
        echo $GENE_START_LOCUS

        GENE_END_LOCUS=$( grep 'GRCh38' gene_info.txt | awk -F '\t' '{print $5}')
        echo $GENE_END_LOCUS

  >>>

    output {
        File results = "gene_info.txt"
        String ASSEMBLY
        String CHR_ID
        Int GENE_START_LOCUS
        Int GENE_END_LOCUS
  }
    runtime {
        memory: memSizeGB + " GB"
        docker: "ubuntu:18.04"
        preemptible: 1
    }

}


