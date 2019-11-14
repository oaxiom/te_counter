'''

Package for genome annotations

'''

from .. import common

valid_assemblies = {'mm10', 'hg38'}
valid_modes = {'genes_tes', 'enhancers'}

class indices:
    def __init__():
        """
        **Purpose**
            Genome container

        """

def build_index(assembly, mode='genes_tes'):
    """
    **Purpose**
        Build an index for te_count.py

    **Arguments**
        assembly (Required)
            build for the specified assembly.

            valid assemblies are {0}

        mode (Optional, default='genes_tes')
            build the index for the specified mode

            valid modes are {1}

    **Returns**
        None

    """.format(valid_assemblies, valid_modes)
