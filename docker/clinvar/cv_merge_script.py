import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Merge clinvar files.')
parser.add_argument('-f', help='basic clinvar base file', required=True)
parser.add_argument('-s', help='trait set file', required=True)
parser.add_argument('-m', help='trait mapping', required=True)
parser.add_argument('-o', default='clinvar_variants.csv',help='output CSV with clinvar variants')


args = parser.parse_args()


Gene_CV_basic = pd.read_csv(args.f, delimiter="\t", 
                            names =["row_type", "VCV_ID", "ClinVar_variant_ID", "variant_type","submissions", "version", 
                                "assembly", "Chr", "start", "stop", "ref", "alt","variant_length", "variant_effect","txpt_hgvsc", 
                                "date_created", "date_updated",  "date_submitted", 
                                "review_status","germline_classification",
                                "onco_review_status","oncogenicity_classification",
                                "som_review_status","somatic_classification",
                                "comment",
                                "functional_category", "FA_comment",
                                "CA_ID", "functional_result"] , header=None, keep_default_na=False)

trait_set = pd.read_csv(args.s, delimiter="\t", 
                       names = ["VCV_ID", "TraitSet_ID","TS_Type" , "Trait_ID","Trait_Type",  "ContributesToAggregateClassification","MG_ID" ] )
trait_set['TraitSet_ID'] = trait_set['TraitSet_ID'].fillna("none")
trait_set[['TraitSet_ID' ]] = trait_set[['TraitSet_ID' ]].astype('str')
trait_set[['Trait_ID' ]] = trait_set[['Trait_ID' ]].astype('str')
trait_set['MG_ID'] = trait_set['MG_ID'].fillna("none")
trait_set['Evidence'] = trait_set['Evidence'].fillna("none")
trait_set['Trait_Type'] = trait_set['Trait_Type'].fillna("none")
trait_set['TS_Type'] = trait_set['TS_Type'].fillna("none")
trait_set=trait_set.sort_values([ "VCV_ID", "MG_ID"])
trait_set['MG_ID'] = trait_set[["VCV_ID", "TraitSet_ID","TS_Type" , "Trait_ID","Trait_Type",  "ContributesToAggregateClassification","MG_ID"]].groupby(["VCV_ID", "TraitSet_ID" ])['MG_ID'].transform(lambda x: '_'.join(x))
trait_set['Trait_ID'] = trait_set[["VCV_ID", "TraitSet_ID","TS_Type" , "Trait_ID","Trait_Type",  "ContributesToAggregateClassification","MG_ID"]].groupby(["VCV_ID", "TraitSet_ID" ])['Trait_ID'].transform(lambda x: '&'.join(x))
trait_set =  trait_set[["VCV_ID", "TraitSet_ID","TS_Type" , "Trait_ID","Trait_Type",  "ContributesToAggregateClassification","MG_ID"]].drop_duplicates()


trait_map = pd.read_csv(args.m, delimiter="\t", 
                       names = ["label", "VCV_ID", "CA_ID", "Trait_Type", "MG_ID", "MG_disease_name"] )
trait_map['MG_ID'] = trait_map['MG_ID'].fillna("none")
trait_map['Trait_Type'] = trait_map['Trait_Type'].fillna("none")
trait_map['MG_disease_name'] = trait_map['MG_disease_name'].fillna("none")
trait_map[['CA_ID' ]] = trait_map[['CA_ID' ]].astype('str')
trait_map = trait_map.sort_values([ "VCV_ID", "MG_ID"])
trait_map['MG_ID'] = trait_map[["label", "VCV_ID", "CA_ID", "Trait_Type", "MG_ID", "MG_disease_name"]].groupby(["label", "VCV_ID","Trait_Type", "CA_ID"])['MG_ID'].transform(lambda x: '_'.join(x))
trait_map['MG_disease_name'] = trait_map[["label", "VCV_ID", "CA_ID", "Trait_Type", "MG_ID", "MG_disease_name"]].groupby(["label", "VCV_ID","Trait_Type", "CA_ID","MG_ID"])['MG_disease_name'].transform(lambda x: '&'.join(x))
trait_map = trait_map[["label", "VCV_ID", "CA_ID", "Trait_Type", "MG_ID", "MG_disease_name"]].drop_duplicates()


trait_comb = pd.merge(trait_set, trait_map, how='outer', on=["VCV_ID", "MG_ID"])
trait_comb['MG_ID'] = trait_comb['MG_ID'].fillna("none")
trait_comb['CA_ID'] = trait_comb['CA_ID'].fillna("none")
trait_comb['MG_disease_name'] = trait_comb['MG_disease_name'].fillna("none")
trait_comb['TraitSet_ID'] = trait_comb['TraitSet_ID'].fillna("none")
trait_comb['Evidence'] = trait_comb['Evidence'].fillna("none")
trait_comb['Trait_Type_y'] = trait_comb['Trait_Type_y'].fillna("none")
trait_comb['Trait_Type_x'] = trait_comb['Trait_Type_x'].fillna("none")
trait_comb['TS_Type'] = trait_comb['TS_Type'].fillna("none")
trait_comb = trait_comb[trait_comb.CA_ID != "none"]
trait_comb = trait_comb.drop(trait_comb[(trait_comb['TraitSet_ID'] == "none") & (trait_comb['Trait_Type_y'] == "Finding")].index)

clinvar_complete = pd.merge(Gene_CV_basic, trait_comb, how='outer', on=["VCV_ID", "CA_ID"])
clinvar_complete['Chr'] = 'chr' + clinvar_complete['Chr'].astype(str)
clinvar_complete['CERFAC_variant_id'] = clinvar_complete[['assembly', 'Chr','start','stop','ref','alt' ]].astype(str).agg(':'.join, axis=1)

cols = ['row_type',
 'VCV_ID','CERFAC_variant_id','ClinVar_variant_ID','submissions',
 #'assembly','Chr','start','stop','ref','alt',
 'variant_length', 'MG_disease_name','ContributesToAggregateClassification',
 'variant_type', 'variant_effect','txpt_hgvsc', 'review_status',
 'germline_classification', 'oncogenicity_classification',
 'somatic_classification',
 'comment',
 'functional_category',
 'FA_comment',
 'functional_result',
 'date_created',
 'date_updated',
 'date_submitted',
'version']

clinvar_complete = clinvar_complete[cols]
clinvar_complete.to_csv(args.o, sep=',', index=False )

