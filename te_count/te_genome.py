#!/usr/bin/env python3

'''

A very simple counter for

'''

import sys, os, argparse, logging
from . import common

# Genome builders:
from .genome import make_genes_tes, make_enh

class index:
    def __init__(self):
        pass

    def build_index(self, genome, mode):
        """
        **Purpose**
            Build an index for te_count

        **Arguments**
            assembly (Required)
                build for the specified assembly.

                valid assemblies are {0}

            mode (Optional, default='genes_tes')
                build the index for the specified mode

                valid modes are {1}

        **Returns**
            None

        """.format(common.valid_assemblies, common.valid_modes)
        assert genome in common.valid_assemblies, '{0} genome assembly not in the list of valid assemblies: {1}'.format(genome, common.valid_assemblies)
        assert mode in common.valid_modes, '{0} mode not in the list of valid modes: {1}'.format(mode, common.valid_modes)

        if mode == 'genes_tes':
            return make_genes_tes(genome)
        elif mode == 'enhancers':
            return make_enh(genome)
        elif more == 'tes': # TODO
            pass

        return False

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt\n")
        sys.exit(0)
