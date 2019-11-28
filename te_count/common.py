
import sys

valid_assemblies = {'mm10', 'hg38'}
valid_modes = {'genes_tes', 'enhancers', 'custom'}

def print_genomes(log=None):
    log.info('  Valid Genome assemblies are:')
    log.info('    hg38 - human')
    log.info('    mm10 - mouse')

