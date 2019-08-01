
import sys

def print_species(log=None):
    log.info('  Valid Species codes are:')
    log.info('    hg38 - human')
    log.info('    mm10 - mouse')

def check_species(species):
    if species not in ('mm10', 'hg38'):
        return False
    return True

