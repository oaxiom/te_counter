'''
Build the combined gencode and repeat data for te_hic

'''

import sys, os

gtf_format = {"feature_type": 1, "feature": 2, "gtf_decorators": 8, "commentlines": "#",
        "loc": "location(chr=column[0], left=column[3], right=column[4])",
        "strand": 6, "skiplines": -1, "force_tsv": True}

rmsk_track_form = {"force_tsv": True, 'loc': 'location(chr=column[5], left=column[6], right=column[7])',
    'strand': 9, 'repName': 10, 'repClass': 11, 'repFamily': 12}

def make_genes_tes(genome, log):
    from .. import miniglbase

    # We have to hardcode the gencode URL as it's location can change and the naming is irregular
    if genome == 'mm10':
        gencode_name = 'gencode.vM23.annotation.gtf.gz'
        gencode_url = 'http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/{0}'.format(gencode_name)
        gene_biotype_label = 'gene_type'
        transcript_biotype_label = 'transcript_type'
    elif genome == 'hg38':
        gencode_name = 'gencode.v42.annotation.gtf.gz'
        gencode_url = 'http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/{0}'.format(gencode_name)
        gene_biotype_label = 'gene_type'
        transcript_biotype_label = 'transcript_type'
    elif genome == 'macFas5':
        gencode_name = 'Macaca_fascicularis.Macaca_fascicularis_6.0.115.gtf.gz'
        gencode_url = 'https://ftp.ensembl.org/pub/release-115/gtf/macaca_fascicularis/{0}'.format(gencode_name)
        gene_biotype_label = 'gene_biotype'
        transcript_biotype_label = 'transcript_biotype'

    repeat_name = f'{genome}_rmsk.txt.gz'
    final_name = f'{genome}_genes_tes.glb'

    script_path = os.path.dirname(os.path.realpath(__file__))

    # Do the downloads, rely on wget to confirm if they are valid downloads:
    os.system(f'wget -c -O {script_path}/{gencode_name}         {gencode_url}')
    os.system(f'wget -c -O {script_path}/{repeat_name}          http://hgdownload.soe.ucsc.edu/goldenPath/{genome}/database/rmsk.txt.gz')
    os.system(f'wget -c -O {script_path}/{genome}.chromSizes.gz http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/database/chromInfo.txt.gz')
    print('Decompressing/Filtering')
    if sys.platform == 'darwin':
        os.system(f"gunzip -c {script_path}/{genome}.chromSizes.gz | grep -v -E 'random|chrUn|chrM|fix|alt'  >{script_path}/{genome}.chromSizes.clean")
    else:
        os.system(f"gunzip -c {script_path}/{genome}.chromSizes.gz | grep -v -E 'random|chrUn|chrM|fix|alt'  >{script_path}/{genome}.chromSizes.clean")

    chr_set = frozenset(['X', 'Y', 'M'] + ['%s' % i for i in range(1, 30)])

    repeats = miniglbase.delayedlist(filename='{0}/{1}'.format(script_path, repeat_name), gzip=True, format=rmsk_track_form)
    gencode = miniglbase.delayedlist('{0}/{1}'.format(script_path, gencode_name), gzip=True, format=gtf_format)

    keep_classes = frozenset(['LINE', 'LTR', 'SINE', 'DNA', 'Retroposon', 'tRNA']) # Retroposon is for human, but we can safely keep it here

    added = 0

    newl = []

    print('Adding Repeats')
    p = miniglbase.progressbar(len(repeats))
    for idx, item in enumerate(repeats):
        if item['repClass'] not in keep_classes:
            continue

        #print(item)
        if str(item['loc']['chr']) not in chr_set:
            continue

        te_name = f"{item['repClass']}:{item['repFamily']}:{item['repName']}"

        newentry = {'loc': item['loc'],
            'strand': item['strand'],
            'name': te_name,
            'type': 'TE',
            'ensg': te_name,
            }
        newl.append(newentry)
        #print(newentry)
        added += 1

        p.update(idx)
    print(f'\nAdded {added:,} features')

    print('Adding Gencode exons')
    added = 0
    p = miniglbase.progressbar(len(gencode))
    for idx, item in enumerate(gencode):
        if item['feature'] != 'exon': # i.e. only include in the annotation if it is an exon
            continue

        if item[gene_biotype_label] not in ('protein_coding', 'lncRNA', 'lincRNA'):
            continue

        if item[transcript_biotype_label] not in ('protein_coding', 'lncRNA', 'lincRNA'):
            continue

        if item['loc']['chr'] not in chr_set:
            continue

        if 'gene_name' in item: # For macFas5 genome;
            gene_name = item['gene_name']
        else:
            gene_name = item['gene_id']

        newentry = {'loc': item['loc'],
            'strand': item['strand'],
            'name': gene_name,
            'type': item[gene_biotype_label],
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

def make_enh(genome, log):
    """
    **Purpose**
        make a genome from the enhancers

    """

    from .. import miniglbase
    script_path = os.path.dirname(os.path.realpath(__file__))

    filename = 'enh_{0}.bed.gz'.format(genome)
    final_name = '{0}_enhancers.glb'.format(genome)

    enh_url = 'http://fantom.gsc.riken.jp/5/datafiles/reprocessed/{0}_latest/extra/enhancer/F5.{0}.enhancers.bed.gz'.format(genome)
    os.system('wget -c -O {0}/{1} {2}'.format(script_path, filename, enh_url))

    script_path = os.path.dirname(os.path.realpath(__file__))

    chr_set = frozenset(['X', 'Y'] + ['%s' % i for i in range(1, 30)])

    enh = miniglbase.genelist(filename=os.path.join(script_path, filename), gzip=True, format={'force_tsv': True, 'loc': 'location(chr=column[0], left=column[1], right=column[2])'})

    newl = []
    added = 0
    p = miniglbase.progressbar(len(enh))
    for idx, item in enumerate(enh):
        newentry = {'loc': item['loc'],
            'name': 'F5enh_{0}_{1}_{2}_{3}'.format(genome, item['loc']['chr'], item['loc']['left'], item['loc']['right']),
            'type': 'enhancer',
            'ensg': 'F5enh_{0}_{1}_{2}_{3}'.format(genome, item['loc']['chr'], item['loc']['left'], item['loc']['right']),
            }

        if item['loc']['chr'] not in chr_set:
            continue

        newl.append(newentry)
        added += 1

    print('\nAdded {:,} features'.format(added))
    p.update(idx)
    gl = miniglbase.genelist()
    gl.load_list(newl)
    gl.save('{0}/{1}'.format(script_path, final_name))

def make_custom(gtffilename, filename_for_index, log):
    """
    **Purpose*
        make a custom genome assembly from gtffilename

    """
    from .. import miniglbase

    chr_set = frozenset(['X', 'Y'] + ['%s' % i for i in range(1, 30)])

    if '.gz' in gtffilename:
        gtf = miniglbase.delayedlist(gtffilename, gzip=True, format=gtf_format)
    else:
        gtf = miniglbase.delayedlist(gtffilename, gzip=False, format=gtf_format)

    # check that we have a 'ensg' key that we can use
    log.info('Checking GTF format')
    for idx, item in enumerate(gtf):
        assert 'ensg' in item, 'the custom GTF must have an "ensg" key that is unique and can be used to score expression'
        assert 'gene_type' in item, 'the custom GTF must have a "gene_type" key, which should be one of "protein_coding/lncRNA"'
        assert 'gene_id' in item, 'the custom GTF must have a "gene_id" key, which should contain a unqie identifier for the gene'
        if idx > 10: break # check the first 10 entries
    gtf.reset()

    added = 0
    newl = []
    log.info('Processing custom gtf {0}'.format(gtffilename))
    p = miniglbase.progressbar(len(gtf))
    for idx, item in enumerate(gtf):
        if item['feature'] != 'exon': # i.e. only include in the annotation if it is an exon
            continue

        if item['loc']['chr'] not in chr_set:
            continue

        # need to deal with GTF keys intelligently...
        newentry = {'loc': item['loc'],
            'strand': item['strand'],
            'name': item['gene_name'],
            'type': item['gene_type'],
            'ensg': item['gene_id'],
            }
        newl.append(newentry)
        added += 1

        p.update(idx)

    print('\nAdded {:,} features'.format(added))

    gl = miniglbase.genelist()
    gl.load_list(newl)
    gl.save(filename_for_index)

    return True

def make_snrnps(genome, log):
    from .. import miniglbase

    # We have to hardcode the gencode URL as its location can change and the naming is irregular
    if genome == 'mm10':
        gencode_name = 'gencode.vM23.annotation.gtf.gz'
        gencode_url = 'http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/{0}'.format(gencode_name)
    elif genome == 'hg38':
        gencode_name = 'gencode.v42.annotation.gtf.gz'
        gencode_url = 'http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/{0}'.format(gencode_name)

    final_name = f'{genome}_snrnps.glb'

    script_path = os.path.dirname(os.path.realpath(__file__))

    # Do the downloads, rely on wget to confirm if they are valid downloads:
    os.system(f'wget -c -O {script_path}/{gencode_name}         {gencode_url}')
    os.system(f'wget -c -O {script_path}/{genome}.chromSizes.gz http://hgdownload.cse.ucsc.edu/goldenPath/{genome}/database/chromInfo.txt.gz')
    print('Decompressing/Filtering')
    if sys.platform == 'darwin':
        os.system(f"gunzip -c {script_path}/{genome}.chromSizes.gz | grep -v -E 'random|chrUn|chrM|fix|alt'  >{script_path}/{genome}.chromSizes.clean")
    else:
        os.system(f"gunzip -c {script_path}/{genome}.chromSizes.gz | grep -v -E 'random|chrUn|chrM|fix|alt'  >{script_path}/{genome}.chromSizes.clean")

    chr_set = frozenset(['X', 'Y'] + ['%s' % i for i in range(1, 30)])

    gencode = miniglbase.delayedlist(f'{script_path}/{gencode_name}', gzip=True, format=gtf_format)

    print('Adding snRNAs')
    newl = []
    p = miniglbase.progressbar(len(gencode))
    for idx, item in enumerate(gencode):
        if item['feature'] != 'exon': # i.e. only include in the annotation if it is an exon
            continue

        if item['gene_type'] not in ('snRNA'):
            continue

        if item['transcript_type'] not in ('snRNA'):
            continue

        if item['loc']['chr'] not in chr_set:
            continue

        newentry = {'loc': item['loc'],
            'strand': item['strand'],
            'name': 'snRNA-' + item['gene_name'],
            'type': item['gene_type'],
            'ensg': item['gene_id'].split('.')[0],
            }
        newl.append(newentry)

        p.update(idx)

    print(f'\nAdded {len(newl):,} features')

    gl = miniglbase.genelist()
    gl.load_list(newl)
    gl.save(f'{script_path}/{final_name}')

    return True
