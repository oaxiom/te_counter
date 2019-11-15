'''

Package for genome annotations

'''

import os

valid_assemblies = {'mm10', 'hg38'}
valid_modes = {'genes_tes', 'enhancers'}

from .. import common
from .make import make_genes_tes

def check_genome_available(genome, mode):
    """
    **Purpose**
        Check that the genome has been built

    **Returns**
        True of False
    """
    assert genome in common.valid_assemblies, '{0} genome assembly not in the list of valid assemblies: {1}'.format(genome, common.valid_assemblies)
    assert mode in common.valid_modes, '{0} mode not in the list of valid modes: {1}'.format(mode, common.valid_modes)

    script_path = os.path.dirname(os.path.realpath(__file__))
    filename = '{0}_{1}.glb'.format(genome, mode)

    # Check the filename is here:
    if os.path.exists(os.path.join(script_path, filename)):
        return True
    return False
