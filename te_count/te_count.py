'''

A counter for various genome features in a range of data.

'''

import sys, os, argparse, logging
from collections import defaultdict
from operator import itemgetter
from . import miniglbase # miniglbase namespace mangling!
import pysam
from . import common

class measureTE:
    def __init__(self, base_path, quality_threshold):
        '''
        **Purpose**
            Constructor

        **Arguments**
            base_path (Required)
                the path we are being run in


        '''
        self.base_path = base_path
        self.total_reads = 0
        self.quality_threshold = quality_threshold

    def bind_genome(self, genelist_glb_filename):
        self.genome = miniglbase.glload(genelist_glb_filename)
        self.all_feature_names = sorted(list(set(self.genome['ensg'])))

    def parse_bampe(self, filename, strand=False, log=None):
        '''
        **Purpose**
            Load in a BAMPE file

        **Arguments**
            filename (Required)
                filename of the BAMPE file
        '''
        assert filename, 'You must specify a filename'

        if strand:
            raise NotImplementedError()

        final_results = {i: 0 for i in self.all_feature_names}

        bucket_size = miniglbase.config.bucket_size

        sam = pysam.AlignmentFile(filename, 'r')
        idx = 0

        try:
            while 1:
                idx += 1
                read1 = next(sam)
                read2 = next(sam)

                if int(read1.mapping_quality) < self.quality_threshold: # these guys share mapq
                    continue

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
                    loc1 = miniglbase.location(chr=chrom, left=left, right=left+1)
                    loc2 = miniglbase.location(chr=chrom, left=rite-1, right=rite) # I just align the two edges, then I don't need to worry about split-reads, and I rely on the duplicate removal
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
                        if 'protein_coding' in types or 'lincRNA' in types or 'lncRNA' in types:
                            for e in ensgs:
                                if ':' in ensgs: # A TE, skip it
                                    continue
                                final_results[e] += 1
                        elif 'TE' in types:
                            for e in ensgs: # Not in any other mRNA, so okay to count as a TE
                                final_results[e] += 1
                        elif 'enhancer' in types:
                            for e in ensgs: # Not in any other mRNA, so okay to count as a TE
                                final_results[e][barcode] += 1
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

    def parse_bamse(self, filename, strand=False, log=None):
        '''
        **Purpose**
            Load in a BAMSE file

        **Arguments**
            filename (Required)
                filename of the BAMSE file
        '''
        assert filename, 'You must specify a filename'

        if strand:
            raise NotImplementedError()

        final_results = {i: 0 for i in self.all_feature_names}

        bucket_size = miniglbase.config.bucket_size

        sam = pysam.AlignmentFile(filename, 'r')
        idx = 0

        try:
            while 1:
                idx += 1
                read1 = next(sam)
                if read1.is_unmapped or read1.is_duplicate or read1.is_qcfail:
                    continue

                if int(read1.mapping_quality) < self.quality_threshold:
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
                    loc1 = miniglbase.location(chr=chrom, left=left, right=left+1)
                    loc2 = miniglbase.location(chr=chrom, left=rite-1, right=rite) # I just align the two edges, then I don't need to worry about split-reads, and I rely on the duplicate removal
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
                        if 'protein_coding' in types or 'lincRNA' in types or 'lncRNA' in types:
                            for e in ensgs:
                                if ':' in ensgs: # A TE, skip it
                                    continue
                                final_results[e] += 1
                        elif 'TE' in types:
                            for e in ensgs: # Not in any other mRNA, so okay to count as a TE
                                final_results[e] += 1
                        elif 'enhancer' in types:
                            for e in ensgs: # Not in any other mRNA, so okay to count as a TE
                                final_results[e][barcode] += 1
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

    def save_result_bulk(self, result, out_filename, log=None):
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

    def sc_parse_bamse(self, filename, UMIS=True, whitelistfileanem=None, strand=False, log=None):
        '''
        **Purpose**
            Load in a BAMSE file, for single cell data, and look for the CR and UMI tags.

        **Arguments**
            filename (Required)
                filename of the BAMSE file

            UMIS (OPtional, default=True)
                Wheter to get the UMIs out of UR tag

            whitelist (Required)
                perform whitlisting on the barcodes;

        '''
        assert filename, 'You must specify a filename'
        assert whitelistfileanem, 'You must specify a filename for the barcode whitelist'

        if not os.path.exists(whitelistfileanem):
            raise AssertionError('{0} -w whitelist file not found'.format(whitelistfileanem))

        whitelist = []
        oh = open(whitelistfileanem)
        for line in oh:
            whitelist.append(line.strip())
        oh.close()
        whitelist = set(whitelist)

        final_results = {i: defaultdict(int) for i in self.all_feature_names} # pseudo-sparse array
        self.barcodes = {}
        umis = defaultdict(set)
        bucket_size = miniglbase.config.bucket_size
        read_assinged_to_gene = 0
        valid_barcodes_reads = 0
        invalid_barcode_reads = 0
        sam = pysam.AlignmentFile(filename, 'r')
        idx = 0

        try:
            while 1:
                idx += 1
                if idx % 1000000 == 0:
                    log.info('Processed {:,} SE reads'.format(idx))
                    #break

                read = next(sam)
                if read.is_unmapped or read.is_duplicate or read.is_qcfail:
                    continue

                if int(read.mapping_quality) < self.quality_threshold:
                    continue

                # Check we have a CR:Z and UR:Z key:
                tags = dict(read.get_tags())
                if 'CR' not in tags: # No barcode
                    continue
                if UMIS and 'UR' not in tags: # No UMI
                    continue

                barcode = tags['CR']
                if barcode not in whitelist:
                    # TODO: 1 bp mismatch recovery
                    invalid_barcode_reads += 1
                    continue

                if UMIS:
                    umi = '{0}-{1}'.format(tags['UR'], barcode) # UMI should be unique for both
                else:
                    umi = barcode # putting this here like this will enforce a policy of 1 read per genomic-location-strand.

                chrom = read.reference_name.replace('chr', '')
                if chrom not in self.genome.buckets: # Must be a valid chromosome
                    continue

                left = read.reference_start
                rite = read.reference_end
                loc_strand = '-' if read.is_reverse else '+'

                # Check we havne't seen this UMI/CB before:
                if umi in umis: # umi/CB was seen
                    l = (chrom, left, right)
                    if UMIS and l in umis[umi]: # check we haven't seen this exact fragment;
                        continue # We've seen this umi and loc before
                    umis[umi].add(l)

                # reach into the genelist guts...
                # work out which of the buckets is required:
                left_buck = ((left-1)//bucket_size) * bucket_size
                right_buck = (rite//bucket_size) * bucket_size
                buckets_reqd = range(left_buck, right_buck+bucket_size, bucket_size)
                result = []
                # get the ids reqd.
                loc_ids = set()
                if buckets_reqd:
                    loc1 = miniglbase.location(chr=chrom, left=left, right=left+1)
                    loc2 = miniglbase.location(chr=chrom, left=rite-1, right=rite) # I just align the two edges, then I don't need to worry about split-reads, and I rely on the duplicate removal
                    # To get rid of the same gene twice
                    for buck in buckets_reqd:
                        if buck in self.genome.buckets[chrom]:
                            loc_ids.update(self.genome.buckets[chrom][buck]) # set = unique ids

                    for index in loc_ids:
                        # check strands

                        if strand:
                            if loc_strand != self.genome.linearData[index]["strand"]:
                                continue

                        #print loc.qcollide(self.linearData[index]["loc"]), loc, self.linearData[index]["loc"]
                        if loc1.qcollide(self.genome.linearData[index]["loc"]): # Any 1 bp overlap...
                            result.append(self.genome.linearData[index])

                        if loc2.qcollide(self.genome.linearData[index]["loc"]): # Any 1 bp overlap...
                            result.append(self.genome.linearData[index])

                    if result:
                        # We are going to add to something:
                        if barcode not in self.barcodes:
                            self.barcodes[barcode] = 0
                        self.barcodes[barcode] += 1
                        # do the annotation so that a read only gets counted to a TE if it does not hit a gene:
                        #for r in result:
                        #    print(r)

                        # This will currently allow 1 read to be counted twice if each edge is inside a different feature.
                        # Is that wrong, or a reasonable compromise?

                        types = set([i['type'] for i in result])
                        ensgs = set([i['ensg'] for i in result]) # only count 1 read to 1 gene
                        if 'protein_coding' in types or 'lincRNA' in types or 'lncRNA' in types:
                            for e in ensgs:
                                if ':' in ensgs: # A TE, skip it
                                    continue
                                final_results[e][barcode] += 1
                        elif 'TE' in types:
                            for e in ensgs: # Not in any other mRNA, so okay to count as a TE
                                final_results[e][barcode] += 1
                        elif 'enhancer' in types:
                            for e in ensgs: # Not in any other mRNA, so okay to count as a enhancer
                                final_results[e][barcode] += 1
                        read_assinged_to_gene += 1
                        #print()

        except StopIteration:
            pass # the last read

        sam.close()
        log.info('Processed {0:,} SE reads'.format(idx))
        log.info('Found {0:,} invalid barcode reads'.format(invalid_barcode_reads))
        log.info('Assigned {0:,} ({1:.1f}%) reads to genes'.format(read_assinged_to_gene, (read_assinged_to_gene/idx * 100.0))) # add per cents here;
        self.total_reads = idx

        return final_results

    def sc_save_result(self, result, out_filename, maxcells, log=None):
        '''
        **Purpose**
            Save the data to a TSV file

        **Arguments**
            out_filename (Required)
                the filename to save the data to
        '''
        assert out_filename, 'You must specify a filename'

        #import matplotlib.pyplot as plot

        log.info('Densifying and saving "{0}"'.format(out_filename))
        log.info('Found {0:,} barcodes'.format(len(self.barcodes)))

        barcodes_to_do = sorted(self.barcodes.items(), key=itemgetter(1), reverse=True)
        if len(self.barcodes) > maxcells: # Or dont bother doing
            # Work out the maxcells barcodes to save
            log.info('Keeping the best {0:,} barcodes'.format(maxcells))
            barcodes_to_do = [i[0] for i in barcodes_to_do][0:maxcells]
        elif maxcells > len(self.barcodes):
            log.warning('Asked for {0:,} maxcells, but only {1:,} barcodes found'.format(maxcells, len(self.barcodes)))
            barcodes_to_do = [i[0] for i in barcodes_to_do]
        else:
            barcodes_to_do = [i[0] for i in barcodes_to_do]

        if '.tsv' not in out_filename:
            out_filename = '{0}.tsv'.format(out_filename)
        barcode_freq_filename = out_filename.replace('.tsv', '.barcode_freq.tsv')

        oh = open(barcode_freq_filename, 'w')
        for b in barcodes_to_do:
            oh.write('{0}\t{1}\n'.format(b, self.barcodes[b]))
        oh.close()

        log.info('Saving barcode read frequency file to {0}'.format(barcode_freq_filename))

        oh = open(out_filename, 'w')
        oh.write('{0}\t{1}\n'.format('name', '\t'.join(result.keys())))

        for barcode in barcodes_to_do:
            counts = []
            for feature in result:
                if barcode in result[feature]: # Stop defaultdict from densifying
                    counts.append(result[feature][barcode])
                else:
                    counts.append(0)
            oh.write('{0}\n'.format('\t'.join([barcode] + [str(c) for c in counts])))
        oh.close()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt\n")
        sys.exit(0)
