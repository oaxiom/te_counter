'''
Build the combined gencode and repeat data for te_hic

'''

import sys, os

def make_genes_tes(genome):
    from .. import miniglbase

    # We have to hardcode the gencode URL as it's location can change and the naming is irregular
    if genome == 'mm10':
        gencode_name = 'gencode.vM23.annotation.gtf.gz'
        gencode_url = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/{0}'.format(gencode_name)
    elif genome == 'hg38':
        gencode_name = 'gencode.v32.annotation.gtf.gz'
        gencode_url = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/{0}'.format(gencode_name)

    repeat_name = '{0}_rmsk.txt.gz'.format(genome)
    final_name = '{0}_genes_tes.glb'.format(genome)

    gtf_format = {"feature_type": 1, "feature": 2, "gtf_decorators": 8, "commentlines": "#",
            "loc": "location(chr=column[0], left=column[3], right=column[4])",
            "strand": 6, "skiplines": -1, "force_tsv": True}

    rmsk_track_form = {"force_tsv": True, 'loc': 'location(chr=column[5], left=column[6], right=column[7])',
        'repName': 10, 'repClass': 11, 'repFamily': 12}

    script_path = os.path.dirname(os.path.realpath(__file__))

    # Do the downloads, rely on wget to confirm if they are valid downloads:
    os.system('wget -c -O {0}/{1}               {2}'.format(script_path, gencode_name, gencode_url))
    os.system('wget -c -O {0}/{1}               http://hgdownload.soe.ucsc.edu/goldenPath/{2}/database/rmsk.txt.gz'.format(script_path, repeat_name, genome))
    os.system('wget -c -O {0}/{1}.chromSizes.gz ftp://hgdownload.cse.ucsc.edu/goldenPath/{1}/database/chromInfo.txt.gz'.format(script_path, genome))
    if sys.platform == 'darwin':
        os.system("gunzip -c {0}/{1}.chromSizes.gz | grep -v -E 'random|chrUn|chrM'  >{0}/{1}.chromSizes.clean".format(script_path, genome))
    else:
        os.system("gunzip -c {0}/{1}.chromSizes.gz | grep -v -r 'random|chrUn|chrM'  >{0}/{1}.chromSizes.clean".format(script_path, genome))

    chr_set = frozenset(['X', 'Y'] + ['%s' % i for i in range(1, 30)])

    repeats = miniglbase.delayedlist(filename='{0}/{1}'.format(script_path, repeat_name), gzip=True, format=rmsk_track_form)
    gencode = miniglbase.delayedlist('{0}/{1}'.format(script_path, gencode_name), gzip=True, format=gtf_format)

    keep_classes = frozenset(['LINE', 'LTR', 'SINE', 'DNA', 'Retroposon']) # Retroposon is for human, but we can safely keep it here

    added = 0

    newl = []

    print('Repeats')
    p = miniglbase.progressbar(len(repeats))
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

    print('\nAdded {:,} features'.format(added))
    print('Gencode')
    p = miniglbase.progressbar(len(gencode))
    for idx, item in enumerate(gencode):
        if item['feature'] != 'exon': # i.e. only include in the annotation if it is an exon
            continue

        if item['gene_type'] not in ('protein_coding', 'lncRNA'):
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

    print('\nAdded {:,} features'.format(added))

    gl = miniglbase.genelist()
    gl.load_list(newl)
    gl.save('{0}/{1}'.format(script_path, final_name))

    return True
