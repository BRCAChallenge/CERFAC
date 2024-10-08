
        import os
        import io
        import json
        import string
        import pandas as pd
        import re
        import logging

        import gnomad
        import hail as hl
        
        hl.init(spark_conf={'spark.driver.memory': "~{hailMemSizeGB}g"})


        
        start_pos =hl.int32("~{GENE_START_LOCUS}")
        stop_pos =hl.int32("~{GENE_END_LOCUS}")

        from gnomad.resources.grch38.gnomad import public_release
        v4exomes = public_release("exomes").ht()
        v4exomes.count()


        from gnomad.resources.grch38.gnomad import public_release
        v4genomes = public_release("genomes").ht()
        v4genomes.count()

        exglobals = v4exomes.select_globals('tool_versions', 'vrs_versions',  'version'  ).head(2)
        flattened = exglobals.globals.flatten()
        flattened.export('gnomad_version_info.tsv')
        badcols = ['variant_effect']
        gnomad_union_df = gnomad_union_df.explode(badcols)


        clinvar_complete = clinvar_complete[clinvar_complete.variant_effect != "splice donor variant"]
        clinvar_complete = clinvar_complete[clinvar_complete.variant_effect != "splice acceptor variant"]
        clinvar_complete = clinvar_complete[clinvar_complete.variant_effect != "intron variant"]


        gnomad_union_df = gnomad_union_df[gnomad_union_df.variant_effect != "splice_region_variant"]
        gnomad_union_df = gnomad_union_df[gnomad_union_df.variant_effect != "splice_donor_region_variant"]
        gnomad_union_df = gnomad_union_df[gnomad_union_df.variant_effect != "splice_donor_variant"]
        gnomad_union_df = gnomad_union_df[gnomad_union_df.variant_effect != "splice_acceptor_variant"]
        gnomad_union_df = gnomad_union_df[gnomad_union_df.variant_effect != "splice_donor_5th_base_variant"]
        gnomad_union_df = gnomad_union_df[gnomad_union_df.variant_effect != "splice_polypyrimidine_tract_variant"]
        gnomad_union_df = gnomad_union_df[gnomad_union_df.variant_effect != "intron_variant"]


#"NM_007294.4(BRCA1):c.5476C&gt;T (p.Gln1826Ter)"
MANE_select_txpt = "NM_007294.4"

√¢‚Ç¨≈°√É‚Äû√É¬∫
√¢‚Ç¨≈°√É‚Äû√É¬π 
√É¬¢√¢‚Äö¬¨√¢‚Ç¨¬π