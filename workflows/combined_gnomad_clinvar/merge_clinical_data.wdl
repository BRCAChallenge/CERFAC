version 1.0

workflow merge_clinical_data {

    meta {
        author: "Allison Cheney"
        email: "archeney@ucsc.edu"
        description: "Add user-supplied clinical data files to gnomad and clinvar variant data"
    }

    parameter_meta {
        GENE_NAME: "Name of relevant gene"
        VARIANTS_FILE: "Gnomad and clinvar variants file"
        FUNCTIONAL_SCORES: "File with functional assay scores and hgvs variant names"
        CLINICAL_DATA: "File with user-supplied case data and hgvs variant names"
    }
    input {
        String GENE_NAME
        File VARIANTS_FILE
        File FUNCTIONAL_SCORES
        File CLINICAL_DATA
    }
    #The order in which the workflow block and task definitions are arranged in the script does not matter. 
    #Nor does the order of the call statements matter, as we'll see further on.

    call merge_variants_clinical{
        input: 
            GENE_NAME=GENE_NAME,
            VARIANTS_FILE=VARIANTS_FILE, 
            CLINICAL_DATA=CLINICAL_DATA, 
            FUNCTIONAL_SCORES=FUNCTIONAL_SCORES
    }


    output{
        File output_clinical_file = merge_variants_clinical.var_clinical_complete    }
}




task merge_variants_clinical {
    input {
        String GENE_NAME
        File VARIANTS_FILE
        File FUNCTIONAL_SCORES
        File CLINICAL_DATA
        Int memSizeGB = 2*round(size(VARIANTS_FILE, "GB") + size(FUNCTIONAL_SCORES, 'GB') + size(CLINICAL_DATA, 'GB')) + 1
        Int threadCount = 1
        Int diskSizeGB = 3*round(size(VARIANTS_FILE, "GB") + size(FUNCTIONAL_SCORES, 'GB') + size(CLINICAL_DATA, 'GB')) + 5

    }

    command <<<
        set -eux -o pipefail
        echo "MEM_SIZE=$MEM_SIZE" >&2
        echo "MEM_UNIT=$MEM_UNIT" >&2

        python3 <<CODE
        import pandas as pd

        variants_df = pd.read_csv("~{VARIANTS_FILE}", delimiter=",", keep_default_na=True)

        def find_file_type_import(file_name):
            read_functions = {'csv': pd.read_csv,
                                        'xlsx': pd.read_excel,
                                        'txt': pd.read_table,
                                        'tsv': pd.read_table}
            [df] = [read(file_name) for file_ext, read in read_functions.items()
                if file_name.endswith(file_ext)]
            return df
        #

        scores_df = find_file_type_import("~{FUNCTIONAL_SCORES}")

        cases_df = find_file_type_import("~{CLINICAL_DATA}")

        clinical_comb_df = pd.merge(scores_df, cases_df, how='outer', on=["hgvs_nt", "hgvs_nt"])

        var_clinical_complete = pd.merge(variants_df, clinical_comb_df, how='outer', on=["hgvs_nt", "hgvs_nt"])

        var_clinical_complete.to_csv("~{GENE_NAME}_variants_clinical.csv", sep=',', index=False )


        CODE

    >>>

    output {
        File var_clinical_complete = "~{GENE_NAME}_variants_clinical.csv"
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








