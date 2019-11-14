#!/usr/bin/env python3

'''

A very simple counter for

'''

import sys, os, argparse, logging
import glbase3 # glbase3 namespace mangling!
import pysam
from . import common

class measureTE:
    def __init__(self, base_path):
        '''
        **Purpose**
            Constructor

        **Arguments**
            base_path (Required)
                the path we are being run in


        '''
        self.base_path = base_path
        self.total_reads = 0

    def bind_genome(self, genelist_glb_filename):
        self.genome = glbase3.glload(genelist_glb_filename)
        self.all_feature_names = sorted(list(set(self.genome['ensg'])))

    def parse_bampe(self, filename, log=None):
        '''
        **Purpose**
            Load in a BAMPE file

        **Arguments**
            filename (Required)
                filename of the BAMPE file
        '''
        assert filename, 'You must specify a filename'

        final_results = {i: 0 for i in self.all_feature_names}

        bucket_size = glbase3.config.bucket_size

        sam = pysam.AlignmentFile(filename, 'r')
        idx = 0

        try:
            while 1:
                idx += 1
                read1 = next(sam)
                read2 = next(sam)
                if read1.is_unmapped or read1.is_duplicate or read1.is_qcfail:
                    continue
                if read2.is_unmapped or read2.is_duplicate or read2.is_qcfail:
                    continue

                if '_'.join(read1.query_name.split('/')[0:-1]) != '_'.join(read2.query_name.split('/')[0:-1]):
                    log.error('Unmatched pair!')
                    sys.quit()

                chrom = read1.reference_name.replace('chr', '')
                left = read1.reference_start
                rite = read2.reference_start

                if chrom not in self.genome.buckets: # Must be a valid chromosome
                    continue

                # reach into the genelist guts...
                # work out which of the buckets is required:
                left_buck = ((left-1)//bucket_size) * bucket_size
                right_buck = (rite//bucket_size) * bucket_size
                buckets_reqd = list(range(left_buck, right_buck+bucket_size, bucket_size))
                result = []
                # get the ids reqd.
                loc_ids = set()
                if buckets_reqd:
                    loc1 = glbase3.location(chr=chrom, left=left, right=left+1)
                    loc2 = glbase3.location(chr=chrom, left=rite-1, right=rite) # I just align the two edges, then I don't need to worry about split-reads, and I rely on the duplicate removal
                    # To get rid of the same gene twice
                    for buck in buckets_reqd:
                        if buck in self.genome.buckets[chrom]:
                            loc_ids.update(self.genome.buckets[chrom][buck]) # set = unique ids

                    for index in loc_ids:
                        #print loc.qcollide(self.linearData[index]["loc"]), loc, self.linearData[index]["loc"]
                        if loc1.qcollide(self.genome.linearData[index]["loc"]): # Any 1 bp overlap...
                            result.append(self.genome.linearData[index])

                        if loc2.qcollide(self.genome.linearData[index]["loc"]): # Any 1 bp overlap...
                            result.append(self.genome.linearData[index])

                    if result:
                        # do the annotation so that a read only gets counted to a TE if it does not hit a gene:
                        #for r in result:
                        #    print(r)
                        types = set([i['type'] for i in result])
                        ensgs = set([i['ensg'] for i in result]) # only count 1 read to 1 gene
                        if 'protein_coding' in types or 'lincRNA' in types:
                            for e in ensgs:
                                if ':' in ensgs: # A TE, skip it
                                    continue
                                final_results[e] += 1
                        elif 'TE' in types:
                            for e in ensgs: # Not in any other mRNA, so okay to count as a TE
                                final_results[e] += 1
                        #print()
                if idx % 1000000 == 0:
                    log.info('Processed {:,} reads'.format(idx))
                    #break
        except StopIteration:
            pass # the last read

        sam.close()
        log.info('Processed {:,} reads'.format(idx))
        self.total_reads = idx

        return final_results

    def parse_bamse(self, filename, log=None):
        '''
        **Purpose**
            Load in a BAMSE file

        **Arguments**
            filename (Required)
                filename of the BAMSE file
        '''
        assert filename, 'You must specify a filename'

        final_results = {i: 0 for i in self.all_feature_names}

        bucket_size = glbase3.config.bucket_size

        sam = pysam.AlignmentFile(filename, 'r')
        idx = 0

        try:
            while 1:
                idx += 1
                read1 = next(sam)
                if read1.is_unmapped or read1.is_duplicate or read1.is_qcfail:
                    continue

                chrom = read1.reference_name.replace('chr', '')
                left = read1.reference_start
                rite = read1.reference_end

                if chrom not in self.genome.buckets: # Must be a valid chromosome
                    continue

                # reach into the genelist guts...
                # work out which of the buckets is required:
                left_buck = ((left-1)//bucket_size) * bucket_size
                right_buck = (rite//bucket_size) * bucket_size
                buckets_reqd = list(range(left_buck, right_buck+bucket_size, bucket_size))
                result = []
                # get the ids reqd.
                loc_ids = set()
                if buckets_reqd:
                    loc1 = glbase3.location(chr=chrom, left=left, right=left+1)
                    loc2 = glbase3.location(chr=chrom, left=rite-1, right=rite) # I just align the two edges, then I don't need to worry about split-reads, and I rely on the duplicate removal
                    # To get rid of the same gene twice
                    for buck in buckets_reqd:
                        if buck in self.genome.buckets[chrom]:
                            loc_ids.update(self.genome.buckets[chrom][buck]) # set = unique ids

                    for index in loc_ids:
                        #print loc.qcollide(self.linearData[index]["loc"]), loc, self.linearData[index]["loc"]
                        if loc1.qcollide(self.genome.linearData[index]["loc"]): # Any 1 bp overlap...
                            result.append(self.genome.linearData[index])

                        if loc2.qcollide(self.genome.linearData[index]["loc"]): # Any 1 bp overlap...
                            result.append(self.genome.linearData[index])

                    if result:
                        # do the annotation so that a read only gets counted to a TE if it does not hit a gene:
                        #for r in result:
                        #    print(r)
                        types = set([i['type'] for i in result])
                        ensgs = set([i['ensg'] for i in result]) # only count 1 read to 1 gene
                        if 'protein_coding' in types or 'lincRNA' in types:
                            for e in ensgs:
                                if ':' in ensgs: # A TE, skip it
                                    continue
                                final_results[e] += 1
                        elif 'TE' in types:
                            for e in ensgs: # Not in any other mRNA, so okay to count as a TE
                                final_results[e] += 1
                        #print()
                if idx % 1000000 == 0:
                    log.info('Processed {:,} SE reads'.format(idx))
                    #break
        except StopIteration:
            pass # the last read

        sam.close()
        log.info('Processed {:,} SE reads'.format(idx))
        self.total_reads = idx

        return final_results

    def save_result(self, result, out_filename, log=None):
        '''
        **Purpose**
            Save the data to a TSV file

        **Arguments**
            out_filename (Required)
                the filename to save the data to
        '''
        assert out_filename, 'You must specify a filename'

        total_reads = self.total_reads/1e6

        oh = open(out_filename, 'w')
        for k in sorted(result.keys()):
            cpm = result[k] / total_reads
            oh.write('{0}\t{1}\t{2}\n'.format(k, result[k], cpm))
        oh.close()
        log.info('Saved {0}'.format(out_filename))

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