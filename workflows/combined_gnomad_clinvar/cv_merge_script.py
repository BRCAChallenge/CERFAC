import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Merge clinvar files.')
parser.add_argument('-f', help='basic clinvar base file', required=True)
parser.add_argument('-s', help='trait set file', required=True)
parser.add_argument('-r', help='list of relevant RCV IDs file', required=True)
parser.add_argument('-c', help='file with RCVs and corresponding CAIDs' , required=True)
parser.add_argument('-t', help='rcv trait file', required=True)
parser.add_argument('-o', default='clinvar_variants.csv',help='output CSV with clinvar variants')


args = parser.parse_args()


Gene_CV_basic = pd.read_csv(args.f, delimiter="\t", 
                            names =["row_type", "VCV_ID", "ClinVar_variant_ID", "variant_type","submissions", "version", 
                                "assembly", "Chr", "start", "stop", "ref", "alt","variant_length", "variant_effect","change",
                                "date_created", "date_updated",  "date_submitted", 
                                "review_status","germline_classification",
                                "onco_review_status","oncogenicity_classification",
                                "som_review_status","somatic_classification","somatic_type", "somatic_significance", "drug_associated",
                                "comment",
                                "functional_category", "FA_comment", "source_type","study_desc",
                                "clinical_assertion_ID",
                                "method", "method_category",
                               "cell_line", "functional_result"] , header=None)
Gene_CV_basic[['version','Chr','start','stop','variant_length','clinical_assertion_ID' ]] = Gene_CV_basic[['version','Chr','start','stop','variant_length','clinical_assertion_ID' ]].astype('Int64')

trait_set = pd.read_csv(args.s, delimiter="\t",  names = ["VCV_ID", "TraitSet_ID", "Trait_ID","ContributesToAggregateClassification"] )

rcv_trait = pd.read_csv(args.t, delimiter="\t", names = ["VCV_ID", "RCV_ID", "TraitSet_ID", "Trait_ID",   "MG_disease_name","MedGen_ID"])

trait_rcv_match = pd.merge(rcv_trait, trait_set, how='outer', on=["VCV_ID", "TraitSet_ID", "Trait_ID" ])

rcv_CAID = pd.read_csv(args.c, delimiter="\t", names = ["VCV_ID","RCV_ID", "clinical_assertion_ID", "SCV_ID"] )

CAID_rcv_match = pd.merge(rcv_CAID, trait_rcv_match, how='outer', on=["VCV_ID", "RCV_ID"])

clinvar_complete = pd.merge(Gene_CV_basic, CAID_rcv_match, how='outer', on=["VCV_ID", "clinical_assertion_ID", "SCV_ID"])
clinvar_complete[['TraitSet_ID','Trait_ID','clinical_assertion_ID' ]] = clinvar_complete[['TraitSet_ID','Trait_ID','clinical_assertion_ID' ]].astype('Int64')


cols = ['row_type',
 'VCV_ID','CERFAC_variant_id','ClinVar_variant_ID','submissions',
 #'assembly','Chr','start','stop','ref','alt',
 'variant_length', 'MG_disease_name','ContributesToAggregateClassification',
 'variant_type', 'variant_effect','change', 'review_status',
 'germline_classification', 'oncogenicity_classification',
 'somatic_classification',
 'somatic_type',
 'somatic_significance',
 'drug_associated',
 'comment',
 'functional_category',
 'FA_comment',
 'source_type',
 'study_desc',
 'method',
 'method_category',
 'cell_line',
 'functional_result',
 'date_created',
 'date_updated',
 'date_submitted',
'version']

clinvar_complete = clinvar_complete[cols]
clinvar_complete.to_csv(args.o, sep=',', index=False )

