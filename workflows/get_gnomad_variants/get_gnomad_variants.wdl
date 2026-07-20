version 1.0

workflow get_gnomad_variants_standalone {

    meta {
        author: "Allison Cheney"
        email: "archeney@ucsc.edu"
        description: "Extract missense gnomad variants for a specified gene (standalone version without ClinVar merge)"
    }

    parameter_meta {
        GENE_NAME: "Name of relevant gene"
    }

    input {
        String GENE_NAME
    }

    call extract_gene_loc {
        input: GENE_NAME=GENE_NAME
    }

    call get_gnomad_variants {
        input:
            GENE_NAME=GENE_NAME,
            CHR_ID=extract_gene_loc.CHR_ID,
            GENE_LENGTH=extract_gene_loc.GENE_LENGTH,
            GENE_START_LOCUS=extract_gene_loc.GENE_START_LOCUS,
            GENE_END_LOCUS=extract_gene_loc.GENE_END_LOCUS
    }

    output {
        File output_gnomad_variants_file = get_gnomad_variants.gnomadvar
        String output_gnomad_variants_count = get_gnomad_variants.gnomad_variants_count
        Int gene_length = extract_gene_loc.GENE_LENGTH
        String chr_id = extract_gene_loc.CHR_ID
        Int gene_start = extract_gene_loc.GENE_START_LOCUS
        Int gene_end = extract_gene_loc.GENE_END_LOCUS
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
        Int GENE_LENGTH = GENE_END_LOCUS - GENE_START_LOCUS
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "brcachallenge/cerfac:clinvar-latest"
        preemptible: 1
    }
}


task get_gnomad_variants {
    input {
        Int memSizeGBbase = 15
        Int threadCount = 1
        Int diskSizeGBbase = 25
        String GENE_NAME
        String CHR_ID
        Int GENE_START_LOCUS
        Int GENE_END_LOCUS
        Int GENE_LENGTH
        File python_script = "get_gnomad_variants.py"
    }

    # Dynamic resource allocation based on gene size
    Int overmilModifier = if GENE_LENGTH >= 1000000 then 30 else 0
    Int overtwomilModifier = if GENE_LENGTH >= 2000000 then 45 else 0

    Int memory_calc = memSizeGBbase + overmilModifier + overtwomilModifier
    Int hailMemSizeGB = floor((0.8* memory_calc) - 2)
    Int diskSizeGB = diskSizeGBbase + overmilModifier + overtwomilModifier + 10

    command <<<
        set -ex -o pipefail
        echo "MEM_SIZE=$MEM_SIZE" >&2
        echo "MEM_UNIT=$MEM_UNIT" >&2

        python3 "~{python_script}" "~{GENE_NAME}" "~{CHR_ID}" ~{GENE_START_LOCUS} ~{GENE_END_LOCUS} ~{hailMemSizeGB}

    >>>

    output {
        File gnomadvar = "~{GENE_NAME}_gnomad_variants_MANE.csv"
        String gnomad_variants_count = read_string("gnomadcount.txt")
    }

    runtime {
        memory: memory_calc + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "brcachallenge/cerfac:gnomad-latest"
        maxRetries: 0
        preemptible: 1
    }
}
