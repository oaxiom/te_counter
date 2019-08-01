
import sys

def print_species():
    print('  Valid Species codes are:')
    print('    hg38 - human')
    print('    mm10 - mouse')

def check_species(species):
    if species not in ('mm10', 'hg38'):
        print('Species "%s" not found' % species)
        print_species()
        print()
        return False
    return True

