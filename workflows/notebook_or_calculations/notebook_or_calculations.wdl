version 1.0

workflow notebook_or_calculations {

    meta {
        description: "Run the OR calculations notebook with papermill on merged VRS variant data"
    }

    parameter_meta {
        merged_variants_file: "Output of merge_clinical_functional_data_vrs_offline (CSV)"
        notebook: "Notebook to execute (notebooks_OR_Calculations_papermill.ipynb)"
        control: "Total control count in the observational data"
        case: "Total case count in the observational data"
        d0: "Population frequency of disease"
        oddsrat: "Odds ratio threshold"
        loss_of_function_threshold: "Functional score at or below which a variant is loss of function (negative values in parentheses in the notebook)"
        functional_threshold: "Functional score at or above which a variant is functional"
    }

    input {
        File merged_variants_file
        File notebook
        Int control = 43481
        Int case = 46306
        Float d0 = 0.12
        Float oddsrat = 5
        Float loss_of_function_threshold = -1.328
        Float functional_threshold = -0.748
    }

    call run_notebook {
        input:
            merged_variants_file = merged_variants_file,
            notebook = notebook,
            control = control,
            case = case,
            d0 = d0,
            oddsrat = oddsrat,
            loss_of_function_threshold = loss_of_function_threshold,
            functional_threshold = functional_threshold
    }

    output {
        File executed_notebook = run_notebook.executed_notebook
    }
}


task run_notebook {
    input {
        File merged_variants_file
        File notebook
        Int control
        Int case
        Float d0
        Float oddsrat
        Float loss_of_function_threshold
        Float functional_threshold
        Int mem_gb = 8
        Int cpu = 2
    }

    command <<<
        set -euo pipefail
        pip install -q papermill
        cp ~{merged_variants_file} input.csv
        papermill ~{notebook} executed_notebook.ipynb \
            -k ir \
            -p input_file input.csv \
            -p control ~{control} \
            -p case ~{case} \
            -p d0 ~{d0} \
            -p oddsrat ~{oddsrat} \
            -p loss_of_function_threshold ~{loss_of_function_threshold} \
            -p functional_threshold ~{functional_threshold}
    >>>

    output {
        File executed_notebook = "executed_notebook.ipynb"
    }

    runtime {
        docker: "quay.io/jupyter/r-notebook@sha256:5f0bac06b517eba722a4775db7a35ff394bbfe11e62d445802a18e942f2c8bcc"
        memory: mem_gb + " GB"
        cpu: cpu
    }
}
