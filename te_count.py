#!/usr/bin/env python3

'''

A very simple counter for

'''

import sys, os
import glbase3 # glbase3 namespace mangling!
import pysam
import common

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

    def bind_genome(self, genelist_glb_filename):
        self.genome = glbase3.glload(genelist_glb_filename)
        print('Loaded %s' % genelist_glb_filename)

    def load_bedpe(self, filename, out_filename):
        '''
        **Purpose**
            Load in a BEDPE file, ideally output by collect_valid_pairs.py, although I guess any valid BEDPE will do

        **Arguments**
            filename (Required)
                filename of the BEDPE file
        '''
        assert filename, 'You must specify a filename'

        done = 0
        bucket_size = glbase3.config.bucket_size

        output = []

        oh = open(filename, 'r')
        for idx, line in enumerate(oh):
            line = line.strip().split('\t')

            # reach into the genelist guts...
            # work out which of the buckets is required:
            loc = glbase3.location(chr=line[0], left=line[1], right=line[2])
            left_buck = int((loc["left"]-1)/bucket_size) * bucket_size
            right_buck = int((loc["right"])/bucket_size) * bucket_size
            buckets_reqd = list(range(left_buck, right_buck+bucket_size, bucket_size))
            result = []
            # get the ids reqd.
            loc_ids = set()
            if buckets_reqd:
                for buck in buckets_reqd:
                    if buck in self.genome.buckets[loc["chr"]]:
                        loc_ids.update(self.genome.buckets[loc["chr"]][buck]) # set = unique ids

                for index in loc_ids:
                    #print loc.qcollide(self.linearData[index]["loc"]), loc, self.linearData[index]["loc"]
                    if loc.qcollide(self.genome.linearData[index]["loc"]):
                        result.append(self.genome.linearData[index])

                read1_feat = []
                read1_type = []
                if result:
                    for r in result:
                        read1_feat.append(r['name'])
                        read1_type.append(r['type'])

            # work out which of the buckets is required:
            loc = glbase3.location(chr=line[3], left=line[4], right=line[5])
            left_buck = int((loc["left"]-1)/bucket_size) * bucket_size
            right_buck = int((loc["right"])/bucket_size) * bucket_size
            buckets_reqd = list(range(left_buck, right_buck+bucket_size, bucket_size))
            result = []
            # get the ids reqd.
            loc_ids = set()
            if buckets_reqd:
                for buck in buckets_reqd:
                    if buck in self.genome.buckets[loc["chr"]]:
                        loc_ids.update(self.genome.buckets[loc["chr"]][buck]) # set = unique ids

                for index in loc_ids:
                    #print loc.qcollide(self.linearData[index]["loc"]), loc, self.linearData[index]["loc"]
                    if loc.qcollide(self.genome.linearData[index]["loc"]):
                        result.append(self.genome.linearData[index])

                read2_feat = []
                read2_type = []
                if result:
                    for r in result:
                        read2_feat.append(r['name'])
                        read2_type.append(r['type'])

            if read1_feat:
                read1_feat = ', '.join(set(read1_feat))
                read1_type = ', '.join(set(read1_type))
            else:
                read1_feat = 'None'
                read1_type = 'None'

            if read2_feat:
                read2_feat = ', '.join(set(read2_feat))
                read2_type = ', '.join(set(read2_type))
            else:
                read2_feat = 'None'
                read2_type = 'None'

            output.append('\t'.join(line[0:3] + [read1_feat, read1_type] + line[3:] + [read2_feat, read2_type]))

            #print(output[-1])

            done += 1

            if done % 1000000 == 0:
                print('Processed: {:,}'.format(done))
                #break

        print('Processed {:,} reads'.format(len(output)))
        oh.close()

        out = open(out_filename, 'w')
        for o in output:
            out.write('%s\n' % o)
        out.close()

    def load_bed(self, filename, out_filename, expand_bed=0):
        '''
        **Purpose**
            Load in a BED file, ideally output by collect_valid_pairs.py,
            although I guess any valid BED will do

            This function is not officially part of te_hic, but could be useful to annotate
            (for example) a BED list of peaks from a ChIP-seq

        **Arguments**
            filename (Required)
                filename of the BED file

            expand_bed (Optional, default=0)
                Optionally expand the BED coordianted left and right by expand_bed

        '''
        assert filename, 'You must specify a filename'

        done = 0
        bucket_size = glbase3.config.bucket_size

        output = []

        oh = open(filename, 'r')
        for idx, line in enumerate(oh):
            line = line.strip().split('\t')

            # reach into the genelist guts...
            # work out which of the buckets is required:
            loc = glbase3.location(chr=line[0], left=int(line[1])-expand_bed, right=int(line[2])+expand_bed)
            left_buck = int((loc["left"]-1)/bucket_size) * bucket_size
            right_buck = int((loc["right"])/bucket_size) * bucket_size
            buckets_reqd = list(range(left_buck, right_buck+bucket_size, bucket_size))
            result = []
            # get the ids reqd.
            loc_ids = set()
            if buckets_reqd:
                for buck in buckets_reqd:
                    if buck in self.genome.buckets[loc["chr"]]:
                        loc_ids.update(self.genome.buckets[loc["chr"]][buck]) # set = unique ids

                for index in loc_ids:
                    #print loc.qcollide(self.linearData[index]["loc"]), loc, self.linearData[index]["loc"]
                    if loc.qcollide(self.genome.linearData[index]["loc"]):
                        result.append(self.genome.linearData[index])

                read1_feat = []
                read1_type = []
                if result:
                    for r in result:
                        read1_feat.append(r['name'])
                        read1_type.append(r['type'])

            if read1_feat:
                read1_feat = ', '.join(set(read1_feat))
                read1_type = ', '.join(set(read1_type))
            else:
                read1_feat = 'None'
                read1_type = 'None'

            output.append('\t'.join(line[0:3] + [read1_feat, read1_type]))

            #print(output[-1])

            done += 1

            if done % 1000000 == 0:
                print('Processed: {:,}'.format(done))
                #break

        print('Processed {:,} reads'.format(len(output)))
        oh.close()

        out = open(out_filename, 'w')
        for o in output:
            out.write('%s\n' % o)
        out.close()

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('\nNot enough arguments')
        print('assign_to_te.py species in.bedpe out.tsv')
        common.print_species()
        print()
        sys.exit()

    species = sys.argv[1]
    if not common.check_species(species):
        sys.exit()

    script_path = os.path.dirname(os.path.realpath(__file__))

    mte = measureTE(sys.argv[0])
    mte.bind_genome(os.path.join(script_path, 'genome/%s_glb_gencode_tes.glb' % species))
    mte.load_bedpe(sys.argv[2], sys.argv[3])


