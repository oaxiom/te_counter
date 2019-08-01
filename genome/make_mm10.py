#!/usr/bin/env python3

'''

Build the combined gencode and repeat data for te_hic

This is for the mm10 genome

'''

from glbase3 import *

rmsk_track_form = {"force_tsv": True, 'loc': 'location(chr=column[5], left=column[6], right=column[7])',
    'repName': 10, 'repClass': 11, 'repFamily': 12}

chr_set = frozenset(['X', 'Y'] + ['%s' % i for i in range(1, 30)])

repeats = delayedlist(filename='mm10_rmsk.txt.gz', gzip=True, format=rmsk_track_form)
gencode = delayedlist('gencode.vM20.annotation.gtf.gz', gzip=True, format=format.gtf)

keep_classes = frozenset(['LINE', 'LTR', 'SINE', 'DNA', 'Retroposon'])

added = 0

newl = []

print('Repeats')
p = progressbar(len(repeats))
for idx, item in enumerate(repeats):
    if item['repClass'] not in keep_classes:
        continue

    #print(item)
    if str(item['loc']['chr']) not in chr_set:
        continue

    newentry = {'loc': item['loc'],
        'name': '%s:%s:%s' % (item['repName'], item['repFamily'], item['repClass']),
        'type': 'TE',
        'ensg': '%s:%s:%s' % (item['repName'], item['repFamily'], item['repClass'])
        }
    newl.append(newentry)
    #print(newentry)
    added += 1

    #if idx > 100000:
    #    break
    p.update(idx)

print('\nAdded %s features' % added)
print('Gencode')
p = progressbar(len(gencode))
for idx, item in enumerate(gencode):
    if item['feature'] != 'exon':
        continue

    #print(item)
    if item['gene_type'] not in ('protein_coding', 'lincRNA'):
        continue

    if item['loc']['chr'] not in chr_set:
        continue

    newentry = {'loc': item['loc'],
        'name': item['gene_name'],
        'type': item['gene_type'],
        'ensg': item['gene_id'].split('.')[0],
        }
    newl.append(newentry)
    added += 1

    #if idx > 100000:
    #    break

    p.update(idx)

print('\nAdded %s features' % added)

gl = genelist()
gl.load_list(newl)
gl.save('mm10_glb_gencode_tes.glb')


