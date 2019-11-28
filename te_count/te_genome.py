#!/usr/bin/env python3

'''

A very simple counter for

'''

import sys, os, argparse, logging
from . import common

# Genome builders:
from .genome import make_genes_tes, make_enh, make_custom

class index:
    def __init__(self):
        pass

    def build_index(self, genome, mode, log, gtffilename=None):
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

            gtffilename (Optional, default=None)
                the filename for the custom gtf when mode=custom

        **Returns**
            None

        """.format(common.valid_assemblies, common.valid_modes)
        assert mode in common.valid_modes, '{0} mode not in the list of valid modes: {1}'.format(mode, common.valid_modes)
        if mode != 'custom':
            assert genome in common.valid_assemblies, '{0} genome assembly not in the list of valid assemblies: {1}'.format(genome, common.valid_assemblies)

        if mode == 'genes_tes':
            return make_genes_tes(genome, log)
        elif mode == 'enhancers':
            return make_enh(genome, log)
        elif mode == 'tes': # TODO
            pass
        elif mode == 'custom':
            return make_custom(gtffilename, genome, log)
        return False

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt\n")
        sys.exit(0)
