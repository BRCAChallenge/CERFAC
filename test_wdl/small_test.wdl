version 1.0





workflow test_workflow {
    meta {
        author: "Allison Cheney"
        email: "archeney@ucsc.edu"
        description: "This is just a test"
    }

    input {
        String GENE_NAME
        String ASSEMBLY
        String CHR_ID
        Int GENE_START_LOCUS
        Int GENE_END_LOCUS
        #Int NS_ENDPOINT = 50
        #Boolean REMOVE_INTRONIC_VAR = true
        #Boolean REMOVE_SPLICING_VAR = true
        #Boolean REMOVE_NS_ENDPOINT_VAR = true
        #Boolean ETHNICITY_2_GROUP = true
        }

    call hello {
        input: 
            GENE_NAME =GENE_NAME,
            ASSEMBLY=ASSEMBLY,
            CHR_ID=CHR_ID,
            GENE_START_LOCUS=GENE_START_LOCUS,
            GENE_END_LOCUS=GENE_END_LOCUS
  }

    output {
        File outputfile = hello.outputfile
        }
}

task hello {
    input {
        String GENE_NAME
        String ASSEMBLY
        String CHR_ID
        Int GENE_START_LOCUS
        Int GENE_END_LOCUS
  }

    command <<<
        echo "hello world" > result_wdl_test.txt
        echo -e  "\nGENE_NAME \tASSEMBLY \tCHR_ID \tGENE_START_LOCUS \tGENE_END_LOCUS" >> result_wdl_test.txt
        echo -e  "\n~{GENE_NAME} \t~{ASSEMBLY} \t~{CHR_ID} \t~{GENE_START_LOCUS} \t~{GENE_END_LOCUS}" >> result_wdl_test.txt

  >>>

    runtime {
        docker: "ubuntu:18.04"
  }

  output {
    File  outputfile = "result_wdl_test.txt"
  }
}
