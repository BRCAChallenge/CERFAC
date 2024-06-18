version 1.0

workflow annotate_functional_variants {

    meta {
        author: "Allison Cheney"
        email: "archeney@ucsc.edu"
        description: "Extract missense gnomad variants for a specified gene"
    }

    parameter_meta {
        GENE_NAME: "Name of relevant gene"
        GNOMAD_SCRIPT:"gnomad_variants_script.py file used for extracting variants from gnomad"
        MERGE_ALL_SCRIPT:"merge_all_variants.py"
    }
    input {
        String GENE_NAME
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
    call merge_variants{
        input:
            MERGE_ALL_SCRIPT=MERGE_ALL_SCRIPT,
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

        echo "getting chromsome of $GENE_NAME"
        awk -F '\t' 'NR == 1 {print $2}' gene_positions.txt | tee CHR_ID
        
        if grep -q -m 1 "GRCh38" gene_positions.txt; then
            ASSEMBLY='GRCh38' #
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
        docker: "allisoncheney/cerfac_terra:clinvar"
        preemptible: 1
    }

}



task get_gnomad_variants {
    input {
        File GNOMAD_SCRIPT
        Int memSizeGB = 4
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




