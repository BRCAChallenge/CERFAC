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
        CLINICAL_SCORES: "File with clinical scores and hgvs variant names"
        CASE_DATA: "File with user-supplied case data and hgvs variant names"
    }
    input {
        String GENE_NAME
        File VARIANTS_FILE
        File CLINICAL_SCORES
        File CASE_DATA
    }
    #The order in which the workflow block and task definitions are arranged in the script does not matter. 
    #Nor does the order of the call statements matter, as we'll see further on.

    call merge_variants_clinical{
        input: 
            GENE_NAME=GENE_NAME,
            VARIANTS_FILE=VARIANTS_FILE, 
            CASE_DATA=CASE_DATA, 
            CLINICAL_SCORES=CLINICAL_SCORES
    }


    output{
        File output_clinical_file = merge_variants_clinical.var_clinical_complete    }
}




task merge_variants_clinical {
    input {
        String GENE_NAME
        File VARIANTS_FILE
        File CLINICAL_SCORES
        File CASE_DATA
        Int memSizeGB = 6
        Int threadCount = 1
        Int diskSizeGB = 3*round(size(VARIANTS_FILE, "GB") + size(CLINICAL_SCORES, 'GB') + size(CASE_DATA, 'GB')) + 5

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

        scores_df = find_file_type_import("~{CLINICAL_SCORES}")

        cases_df = find_file_type_import("~{CASE_DATA}")

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
        import os
        import io
        import json
        import string
        import pandas as pd



        cv_table = pd.read_csv("~{clinvar_var}", sep=',' )
        gnomad_vars = pd.read_csv("~{gnomadvar}", sep=',' )
        gnomad_vars = gnomad_vars.rename(columns={"txpt_hgvsc_gnomad": "hgvs_nt"}, errors='raise')
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

        combined.to_csv( "~{GENE_NAME}_combined_variants.csv",  sep=',', index=True )



        CODE

    >>>

    output {
        File combined_var = "~{GENE_NAME}_combined_variants.csv"
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




