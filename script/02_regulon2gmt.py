import sys
import pandas as pd
from pyscenic.cli.utils import load_signatures

# regulon_file = "../output/01-step2_reg.tsv"
regulon_file = "../output/01-step2_reg.tsv"

project = "../output/02-ifnb_pbmc"

min_regulon_size = 10

def get_motif_logo(regulon):
    base_url = "http://motifcollections.aertslab.org/v10nr_clust/logos/"
    for elem in regulon.context:
        if elem.endswith('.png'):
            return(base_url + elem)

# regulons to gmt file
print('''
##############################################
    1. Transform regulons to gmt file ...
##############################################
    ''')
sys.stdout.flush()
regulons = load_signatures(regulon_file)
select_cols = [i.name for i in regulons if len(i.genes) >= min_regulon_size]
gmt_file = project + ".regulons.gmt"
txt_file = project + ".regulons.txt"
fo1 = open(gmt_file, 'w')
fo2 = open(txt_file, 'w')
for i in regulons:
    if i.name in select_cols:
        motif = get_motif_logo(i)
        genes = "\t".join(i.genes)
        tf = "%s(%sg)" % (i.transcription_factor, len(i.genes))
        fo1.write("%s\t%s\t%s\n" % (tf, motif, genes))
        fo2.write("%s\t%s\t%s\n" % (tf, motif, genes.replace("\t",",")))

