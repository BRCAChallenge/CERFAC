version 1.0

workflow merge_clinical_data {

    meta {
        author: "Allison Cheney"
        email: "archeney@ucsc.edu"
        description: "Add user-supplied clinical data files to gnomad and clinvar variant data"
    }

    parameter_meta {
        GENE_NAME: "Name of relevant gene"
        VARIANTS_FILE: "The gnomad and clinvar calibration variants file generated by previous workflows."
        FUNCTIONAL_SCORES: "File with functional assay scores VCF or HGVS cDNA or HGVS genomic variant names"
        CLINICAL_DATA: "File with user-supplied case data and VCF or HGVS cDNA or HGVS genomic variant names"
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
            FUNCTIONAL_SCORES=FUNCTIONAL_SCORES,
            CLINICAL_DATA=CLINICAL_DATA
    }


    output{
        File output_combined_variants_file = merge_variants_clinical.var_clinical_complete    }
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
        import csv
        import requests
        import re
        from enum import Enum
        import json


        class AnnotationLayer(str, Enum):
            """Create enum for supported annotation layers"""

            HGVS_PROTEIN = "p"
            HGVS_CDNA = "c"
            HGVS_GENOMIC = "g"
        def get_col_name(datafile_name):
            with open (file=datafile_name, mode='rt', newline='') as file:
                rows = find_file_type_import(file, datafile_name)
                headers = next(rows)
                col_name = parse_header(headers)
            return col_name

        def parse_header(header):
            col_name=str(header[0])
            if "\ufeff" in col_name:
                col_name = col_name.replace("\ufeff", "")
            return col_name

        def find_file_type_import(file, file_name):
            read_functions = {'csv': ',', 'xlsx': '\t', 'txt': '\t', 'tsv': '\t'}
            print(file_name) 
            [rows] = [csv.reader(file, delimiter=delim) for file_ext, delim in read_functions.items()
                    if file_name.endswith(file_ext)]
            return rows

        def load_table(datafile_name, url):
            table = {}
            with open (file=datafile_name, mode='rt', newline='') as file:
                rows = find_file_type_import(file, datafile_name)

                csv_headings = next(rows)
                col_name = parse_header(csv_headings)
                print("col_name =", col_name)

                first_row = next(rows)
                first_var= first_row[0]
                first_var=str(first_var.strip())

                if not first_contains_na(first_var):
                    classification = classify_var(first_var)
                    if classification:
                        print(classification)
                        old, new = parse_rows(first_var, url, 1, classification)  
                        store_entry(table,old,new)


                        for n, row in enumerate(rows, start=2):
                            old_coding_id = row[0]
                            if old_coding_id not in table.keys():
                                old, new = parse_rows(old_coding_id, url, n, classification)  
                                store_entry(table,old,new)
                        return table
                    else:
                        print("Error in variant format")
                else:
                    print("Error in variant format.Contains an NA or None")

                    
                    
                


            
        VCF_string = re.compile(
                r"^(chr|chromosome)?(?P<chromosome>([1-9]|[1][0-9]|[2][0-2]|X|Y))-"
                r"(?P<pos>[1-9]\d*)-(?P<ref>[actg]+)-(?P<alt>[actg]+)$",
                re.IGNORECASE,
            )

        under_VCF_string = re.compile(
                r"^(chr|chromosome)?(?P<chromosome>([1-9]|[1][0-9]|[2][0-2]|X|Y))_"
                r"(?P<pos>[1-9]\d*)_(?P<ref>[actg]+)_(?P<alt>[actg]+)$",
                re.IGNORECASE,
            )

        col_VCF_string = re.compile(
                r"^(chr|chromosome)?(?P<chromosome>([1-9]|[1][0-9]|[2][0-2]|X|Y)):"
                r"(?P<pos>[1-9]\d*):(?P<ref>[actg]+):(?P<alt>[actg]+)$",
                re.IGNORECASE,
            )

        HGVS_string = re.compile(
                r"^(?P<accession>(NC_|NM_|NP_|ENSP|ENST)[^:\s]+):(?P<coordinate>[cgp])\.(?P<pos>[1-9]\d*[_]*[1-9]\d*)(?P<change>\S+)$"
            )


        short_HGVS_string = re.compile(
                r"^(?P<coordinate>[cgp])\.(?P<pos>[1-9]\d*)(?P<change>\S+)$"
            )

        simple_HGVS_string = re.compile(
                r"(?P<coordinate>[cgp])\.(?P<pos>[1-9]\d*)(?P<change>\S+)$"
            )

        def classify_var(first_var):
            VCFmatch = VCF_string.match(first_var)
            HGVSmatch = HGVS_string.match(first_var)

            if HGVSmatch:
                match_dict = HGVSmatch.groupdict()

                accession=match_dict["accession"]
                HGVS_type=AnnotationLayer(match_dict["coordinate"])
                change=match_dict["change"]
                classification = str(HGVS_type)

                return classification
            elif VCFmatch:
                classification = "gnomadVCF"
                return classification
            elif first_var.startswith("~{GENE_NAME}"):
                print("Error in ", first_var, "HGVS variants cannot start with gene name" )
                return None
            
            else:
                HGVSmatchsimp = simple_HGVS_string.match(first_var)
                VCF_underscore = under_VCF_string.match(first_var)
                VCF_colon = col_VCF_string.match(first_var)

                if HGVSmatchsimp:
                    missing_colon = first_var.find(":")
                    has_dash = first_var.find("-")
                    HGVSmatchshort = short_HGVS_string.match(first_var)
                    if HGVSmatchshort:
                        print("Error in ", first_var, ". The variant is missing the accession in variant ID" )
                        return None
                    elif missing_colon== -1:
                        print("Error in ", first_var, ". Missing colon in variant ID" )
                        return None
                    elif has_dash != -1:
                        print("Error in ", first_var, ". HGVS variants cannot contain dashes: intronic variants are not supported" )   
                        return None
                    else:
                        print("Error in ", first_var, ". The variant was detected as a HGVS variant but cannot parse it." )
                        return None
                elif VCF_underscore:
                    print("Error in ", first_var, ". VCF variants must contain - not underscores" )
                    return None
                elif VCF_colon:
                    print("Error in ", first_var, ". VCF variants must contain - not colons" )
                    return None
                else:
                    print("no match found for first variant. Are you sure you used a gnomAD VCF or HGVS variant?")
                    return None



        def first_contains_na(old_VCF):
            patterns= ["^None$", r"^-$", "^unknown$" ]
            patternloc = max([old_VCF.find(pattern) for pattern in patterns])
            return True if patternloc >= 0 else False

        def contains_na(old_VCF):
            old_VCF=str(old_VCF)
            patterns= ["None",  "unknown", "NA" ]
            patternloc = max([old_VCF.find(pattern) for pattern in patterns])
            onlydash= re.compile(r"^\-+$")
            varonlydash = re.match(onlydash, old_VCF)
            if patternloc != -1: 
                return True
            elif varonlydash:
                return True
            else:
                return False


        def parse_rows(old_coding_id, url, n, classification):

            try:
                if not contains_na(old_coding_id):
                    vrs_id = get_api_response(url, old_coding_id, n, classification)
                    return old_coding_id, vrs_id
                        
                elif contains_na(old_coding_id):
                    print(n,  old_coding_id,   "is NA, skipped")
                    return old_coding_id, "NA"
                else:
                    print("weird error, skipped?")  
                    
            except IndexError:
                pass


        def store_entry(table, key, value):
            table[key] = value
            
        def var_checker(var, classification):            
            HGVSmatch = HGVS_string.match(var)
            VCFmatch = VCF_string.match(var)

            if classification == "AnnotationLayer.HGVS_GENOMIC":
                if HGVSmatch:
                    match_dict = HGVSmatch.groupdict()

                    accession=match_dict["accession"]
                    ignore=AnnotationLayer(match_dict["coordinate"])
                    change=match_dict["change"]


                    invalid_gpatterns= [re.compile(r"inv"), re.compile(r"^-$"), re.compile(r"\+"), re.compile(r"[BDEFHIJKLMNOPQRSUVWXYZ]+$"), re.compile(r"(del|dup)[ATCG]+")]
                    patternloc = [ re.search(pattern, var) for pattern in invalid_gpatterns]
                    if patternloc:
                        print("The API cannot process HGVS variants with inv, variants with plus or minus signs, \
                            or variants with del or dup followed by bases, or bases other than ATCG")
                    else:
                        print("no invalid pattern")

                        
                else:
                    print("HGVS genomic class was detected but there's something wrong. Maybe invalid character?")
            elif classification == "AnnotationLayer.HGVS_PROTEIN" and HGVSmatch:
                match_dict = HGVSmatch.groupdict()

                accession=match_dict["accession"]
                ignore=AnnotationLayer(match_dict["coordinate"])
                change=match_dict["change"]
                
                print(change)
                    
                

                return change
            elif classification == "AnnotationLayer.HGVS_CDNA":
                if HGVSmatch:
                    match_dict = HGVSmatch.groupdict()

                    accession=match_dict["accession"]
                    ignore=AnnotationLayer(match_dict["coordinate"])
                    change=match_dict["change"]

                    invalid_gpatterns= [re.compile(r"inv"), re.compile(r"dup"), re.compile(r"-"), re.compile(r"\+"),re.compile(r"\*"), re.compile(r"[BDEFHIJKLMNOPQRSUVWXYZ]+$"), re.compile(r"(del|dup)[ATCG]+")]
                    patternloc = [ re.search(pattern, var) for pattern in invalid_gpatterns]
                    if patternloc:
                        print("The API cannot process HGVS coding variants with dup, inv, variants with asterisk, plus or minus signs, or variants with del followed by bases, or variants with bases other than ATCG")
                    else:
                        print("no invalid pattern")
                else:
                    print("HGVS cdna class was detected but there's something wrong. Maybe invalid character?")
            elif classification == "gnomadVCF" and VCFmatch:
                match_dict = VCFmatch.groupdict()
                chromosome = match_dict["chromosome"].upper()
                pos = int(match_dict["pos"])
                ref = match_dict["ref"].upper()
                alt = match_dict["alt"].upper()
                print("Error in variant format. VCF variant detected but API cannot process." )
            
            else:
                print("Error in variant format, cause unknown")


        def get_api_response(url, old_coding_id, n, classification ):
            
            payload = {'q': old_coding_id}
            response = requests.get(url, params=payload)
            if response.status_code == 200:
                data = response.json()
                warning = data["warnings"]
                
                if not warning:
                    new_VRS = data["variations"][0]["id"]
                    print(new_VRS)
                    return new_VRS
                elif warning:
                    print("printing classification", classification)
                    var_checker(old_coding_id, classification)
                    print("Problem with variant:", old_coding_id)
                    print('Error:', warning[0])
                        
                    print("skipping...")
                else:
                    print("Unknown API error")
            else: 
                print(n, old_coding_id)
                print(f"Error: {response.status_code}")
                        
        def merge_table_entries_main(datafile_name, table):
            items = table.items()
            VCF_VRS = pd.DataFrame({'VCF_genomic_ID': [i[0] for i in items], 'VRS_ID': [i[1] for i in items]})
            calib_df = pd.read_csv(datafile_name, delimiter=",", keep_default_na=True, index_col=False, low_memory=False)
            calib_VRS_df = pd.merge(calib_df, VCF_VRS, how='outer', on=["VCF_genomic_ID", "VCF_genomic_ID"])
            return calib_VRS_df
        def find_file_type_import_pandas(file_name):
            read_functions = {'csv': ',', 'xlsx': '\t', 'txt': '\t', 'tsv': '\t'}
            print(file_name) 
            [orig_df] = [pd.read_csv(file_name, sep=delim, keep_default_na=True, index_col=False, low_memory=False) for file_ext, delim in read_functions.items()
                    if file_name.endswith(file_ext)]
            return orig_df
        def merge_table_entries(datafile_name, table):
            col_name= get_col_name(datafile_name)
            items = table.items()
            VCF_VRS = pd.DataFrame({col_name: [i[0] for i in items], 'VRS_ID': [i[1] for i in items]})
                
            orig_df = find_file_type_import_pandas(datafile_name)
            new_VRS_df = pd.merge(orig_df, VCF_VRS, how='outer', on=[col_name, col_name])
            return new_VRS_df

        URL = "https://normalize.cancervariants.org/variation/to_vrs"

        VRS_cases_table = load_table("~{CLINICAL_DATA}", URL)
        Cases_VRS_table = merge_table_entries("~{CLINICAL_DATA}", VRS_cases_table)

        VRS_scores_table = load_table("~{FUNCTIONAL_SCORES}", URL)
        Scores_VRS_table = merge_table_entries("~{FUNCTIONAL_SCORES}", VRS_scores_table)

        
        VRS_var_table = load_table("~{VARIANTS_FILE}", URL)
        Calib_VRS_table = merge_table_entries_main("~{VARIANTS_FILE}", VRS_var_table)


        clinical_comb_df = pd.merge(Scores_VRS_table, Cases_VRS_table, how='outer', on=["VRS_ID", "VRS_ID"])

        var_clinical_complete = pd.merge(variants_df, clinical_comb_df, how='outer', on=["VRS_ID", "VRS_ID"])

        var_clinical_complete.to_csv("~{GENE_NAME}_variants_functional_clinical.csv", sep=',', index=False )


        CODE

    >>>

    output {
        File var_clinical_complete = "~{GENE_NAME}_variants_functional_clinical.csv"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "allisoncheney/cerfac_terra:merge_clinical_data"
        maxRetries: 3
        preemptible: 1
    }
}







