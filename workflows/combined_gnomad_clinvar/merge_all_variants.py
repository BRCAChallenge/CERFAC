import argparse
import os
import io
import json
import string
import pandas as pd


parser = argparse.ArgumentParser(description='Merge gnomad and clinvar variant files.')
parser.add_argument('-g', help='gnomad variants file', required=True)
parser.add_argument('-c', help='clinvar variants file', required=True)
parser.add_argument('-o', default='combined_variants.csv',help='output CSV with both variant types')
args = parser.parse_args()



cv_table = pd.read_csv(args.c, sep=',' )
gnomad_vars = pd.read_csv(args.g, sep=',' )
combined = gnomad_vars.set_index('txpt_hgvsc').join(cv_table.set_index('txpt_hgvsc'), how='outer', lsuffix='_gnomad', rsuffix='_clinvar' )
combined.sort_values(['txpt_hgvsc'])
combined.to_csv( args.o,  sep=',', index=True )
