'''

Package for genome annotations

'''

import os

valid_assemblies = {'mm10', 'hg38', 'macFas5'}
valid_modes = {'genes_tes', 'enhancers', 'snrnps'}

from .. import common
from .make import make_genes_tes, make_enh, make_custom, make_snrnps

def check_genome_available(genome, mode):
    """
    **Purpose**
        Check that the genome has been built

    **Returns**
        True of False
    """

    assert mode in common.valid_modes, f'{mode} mode not in the list of valid modes: {common.valid_modes}'
    if mode != 'custom':
        assert genome in common.valid_assemblies, '{0} genome assembly not in the list of valid assemblies: {1}'.format(genome, common.valid_assemblies)

    script_path = os.path.dirname(os.path.realpath(__file__))
    filename = f'{genome}_{mode}.glb'

    # Check the filename is here:
    return bool(os.path.exists(os.path.join(script_path, filename)))
