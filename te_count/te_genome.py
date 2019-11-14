#!/usr/bin/env python3

'''

A very simple counter for

'''

import sys, os, argparse, logging
import glbase3 # glbase3 namespace mangling!
import pysam
from . import common

class measureTE:
    def bind_genome(self, genelist_glb_filename):
        self.genome = glbase3.glload(genelist_glb_filename)
        self.all_feature_names = sorted(list(set(self.genome['ensg'])))


# Command-line options;
def prepare_parser():
    exmp = 'te_count.py -i in.bam -o out.bam -g genome'

    description = 'Counts up the number of reads that overlap some set of gene/TE features'

    parser = argparse.ArgumentParser(prog='te_count', description=description, epilog=exmp)

    # Optional:
    optional = parser._action_groups.pop()
    optional.add_argument('--se', action='store_true', required=False, help='Set mode to SE (single-end) mode, default is paired-end mode')

    required = parser.add_argument_group('required arguments')

    required.add_argument('-i', '--inbam', nargs=1, required=True, help='the BAM alignment file containing the reads')
    required.add_argument('-o', '--outtsv', nargs=1, required=True, help='the TSV file to save the genes and count data to')
    required.add_argument('-g', '--genome', nargs=1, required=True, help='A txt file to save the observed barcode whitelist to')

    parser._action_groups.append(optional)

    logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)-8s: %(message)s',
                    datefmt='%m-%d %H:%M')

    parser.log = logging.getLogger('te_count')

    return parser

def main():
    assert sys.version_info >= (3, 6), 'Python >=3.6 is required'

    script_path = os.path.dirname(os.path.realpath(__file__))
    parser = prepare_parser()
    args = parser.parse_args()

    log = parser.log

    log.info('Arguments:')
    log.info('  inbam: %s' % args.inbam[0])
    log.info('  outtsv: %s' % args.outtsv[0])
    log.info('  genome: "%s"' % args.genome[0])
    log.info('  single-end mode: {0} (default is PE)'.format(args.se))

    log.info('Finding Genome')
    species = args.genome[0]
    if not common.check_species(species):
        log.error('{0} genome not found'.format(species))
        common.print_species(log=log)
        sys.exit()

    mte = measureTE(sys.argv[0])
    mte.bind_genome(os.path.join(script_path, 'genome/%s_glb_gencode_tes.glb' % species))
    log.info('Found and loaded {0} genome'.format(args.genome[0]))

    if args.se:
        result = mte.parse_bamse(args.inbam[0], log=log)
    else:
        result = mte.parse_bampe(args.inbam[0], log=log)

    mte.save_result(result, args.outtsv[0], log=log)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt\n")
        sys.exit(0)
