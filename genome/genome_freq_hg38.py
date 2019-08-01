#!/usr/bin/env python3
'''

measure the % and number ofbase pairs each type of TE occupuies in the genome

'''

from glbase3 import *

genome = glload('hg38_glb_gencode_tes.glb')

hg38_genome_size = 3096649726 # http://asia.ensembl.org/Homo_sapiens/Info/Annotation

tes = {}

for te in genome:
    if '?' in te:
        continue
    if ':' not in te['name']:
        continue # omit the genes
    if te['name'] not in tes:
        tes[te['name']] = 0
    tes[te['name']] += len(te['loc'])

newl = []
for k in tes:
    newe = {'name': k,
        'genome_count': tes[k],
        'genome_percent': tes[k] / hg38_genome_size * 100.0}
    newl.append(newe)

gl = genelist()
gl.load_list(newl)
gl.sort('name')
gl.saveTSV('hg38_te_genome_freqs.tsv', key_order=['name', 'genome_count', 'genome_percent'])
gl.save('hg38_te_genome_freqs.glb')


