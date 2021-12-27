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

        __read_assinged_to_gene = 0
        __quality_trimmed = 0
        __read_qc_fail = 0
        __invalid_chromosome = 0

        if strand:
            raise NotImplementedError()

        final_results = {i: 0 for i in self.all_feature_names}

        bucket_size = miniglbase.config.bucket_size

        sam = pysam.AlignmentFile(filename, 'r')
        idx = 0

        self_genome_linearData = self.genome.linearData

        # preprocess loc lookups
        loc_lookups = []
        for feature in self_genome_linearData:
            loc_lookups.append((feature['loc']['left'], feature['loc']['right']))

        try:
            while 1:
                idx += 1
                read1 = next(sam)
                read2 = next(sam)

                if read1.is_unmapped or read1.is_duplicate or read1.is_qcfail:
                    __read_qc_fail += 1
                    continue
                if read2.is_unmapped or read2.is_duplicate or read2.is_qcfail:
                    __read_qc_fail += 1
                    continue

                if int(read1.mapping_quality) < self.quality_threshold: # these guys share mapq
                    __quality_trimmed += 1
                    continue

                if '_'.join(read1.query_name.split('/')[0:-1]) != '_'.join(read2.query_name.split('/')[0:-1]):
                    log.error('Unmatched pair!')
                    sys.quit()

                chrom = read1.reference_name.replace('chr', '')
                loc1 = read1.reference_start
                loc2 = read2.reference_start

                if chrom not in self.genome.buckets: # Must be a valid chromosome
                    __invalid_chromosome += 1
                    continue

                # reach into the genelist guts...
                # work out which of the buckets is required:
                left_buck = ((loc1-1)//bucket_size) * bucket_size
                right_buck = ((loc2+1)//bucket_size) * bucket_size
                buckets_reqd = list(set([left_buck, right_buck])) # list(range(left_buck, right_buck+bucket_size, bucket_size))
                result = []

                # I just align the two edges, then I don't need to worry about split-reads, and I rely on the duplicate removal
                # To get rid of the same gene twice
                loc_ids = set()
                for buck in buckets_reqd:
                    if buck in self.genome.buckets[chrom]:
                        loc_ids.update(self.genome.buckets[chrom][buck]) # set = unique ids

                for index in loc_ids:
                    locG_l = loc_lookups[index][0]
                    locG_r = loc_lookups[index][1]

                    if loc1 >= locG_l and loc1+1 <= locG_r:
                        result.append(self_genome_linearData[index])

                    if loc2-1 >= locG_l and loc2 <= locG_r: # Any 1 bp overlap...
                        result.append(self_genome_linearData[index])

                if result:
                    # do the annotation so that a read only gets counted to a TE if it does not hit a gene:
                    #for r in result:
                    #    print(r)
                    types = set([i['type'] for i in result])
                    ensgs = set([i['ensg'] for i in result])
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
                    __read_assinged_to_gene += 1

                if idx % 1e6 == 0:
                    log.info('Processed {:,} reads'.format(idx))

        except StopIteration:
            pass # the last read

        sam.close()
        log.info('Processed {:,} reads'.format(idx))
        log.info('{:,} Reads were assigned to a gene'.format(__read_assinged_to_gene))
        log.info('{:,} Read quality is too low (<{})'.format(__quality_trimmed, self.quality_threshold))
        log.info('{:,} Reads mapped to an invalid chromosome'.format(__invalid_chromosome))
        log.info('{:,} Reads are QC fails'.format(__read_qc_fail))
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

        self_genome_linearData = self.genome.linearData

        # preprocess loc lookups
        loc_lookups = []
        for feature in self_genome_linearData:
            loc_lookups.append((feature['loc']['left'], feature['loc']['right']))

        try:
            while 1:
                idx += 1
                read1 = next(sam)
                if read1.is_unmapped or read1.is_duplicate or read1.is_qcfail:
                    continue

                if int(read1.mapping_quality) < self.quality_threshold:
                    continue

                chrom = read1.reference_name.replace('chr', '')
                loc1 = read1.reference_start
                loc2 = read1.reference_end

                if chrom not in self.genome.buckets: # Must be a valid chromosome
                    continue

                # reach into the genelist guts...
                # work out which of the buckets is required:
                left_buck = ((loc1-1)//bucket_size) * bucket_size
                right_buck = ((loc2+1)//bucket_size) * bucket_size
                buckets_reqd = [left_buck, right_buck] # list(range(left_buck, right_buck+bucket_size, bucket_size))
                result = []

                # get the ids reqd.
                loc_ids = set()
                for buck in buckets_reqd:
                    if buck in self.genome.buckets[chrom]:
                        loc_ids.update(self.genome.buckets[chrom][buck]) # set = unique ids

                for index in loc_ids:
                    locG_l = loc_lookups[index][0]
                    locG_r = loc_lookups[index][1]

                    if loc1 >= locG_l and loc1+1 <= locG_r:
                        result.append(self_genome_linearData[index])

                    if loc2-1 >= locG_l and loc2 <= locG_r: # Any 1 bp overlap...
                        result.append(self_genome_linearData[index])

                if result:
                    # do the annotation so that a read only gets counted to a TE if it does not hit a gene:
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


                if idx % 1e6 == 0:
                    log.info('Processed {:,} SE reads'.format(idx))

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

        with open(out_filename, 'w') as oh:
            for k in sorted(result.keys()):
                cpm = result[k] / total_reads
                oh.write('{0}\t{1}\t{2}\n'.format(k, result[k], cpm))
        log.info('Saved {0}'.format(out_filename))

    def sc_parse_bamse(self, filename, UMIS=True, whitelistfilename=None, strand=False, log=None):
        '''
        **Purpose**
            Load in a BAMSE file, for single cell data, and look for the CR and UMI tags.

        **Arguments**
            filename (Required)
                filename of the BAMSE file

            UMIS (OPtional, default=True)
                Whether to get the UMIs out of UR tag

            whitelist (Required)
                perform whitlisting on the barcodes;

        '''
        assert filename, 'You must specify a filename'

        whitelist = None
        if whitelistfilename:
            if not os.path.exists(whitelistfilename):
                raise AssertionError('{0} -w whitelist file not found'.format(whitelistfilename))
            whitelist = []
            oh = open(whitelistfilename)
            for line in oh:
                whitelist.append(line.strip())
            oh.close()
            whitelist = set(whitelist)

        final_results = {i: {} for i in self.all_feature_names} # pseudo-sparse array
        self.barcodes = {}
        umis = defaultdict(set)
        bucket_size = miniglbase.config.bucket_size
        __read_assinged_to_gene = 0
        valid_barcodes_reads = 0
        __invalid_barcode_reads = 0
        __quality_trimmed = 0
        __read_qc_fail = 0
        __already_seen_umicb = 0
        sam = pysam.AlignmentFile(filename, 'r')
        idx = 0

        self_genome_linearData = self.genome.linearData

        # preprocess loc lookups
        loc_lookups = []
        for feature in self_genome_linearData:
            loc_lookups.append((feature['loc']['left'], feature['loc']['right']))

        try:
            while 1:
                idx += 1
                if idx % 1000000 == 0:
                    log.info('Processed {:,} SE reads'.format(idx))
                    #break

                read = next(sam)
                if read.is_unmapped or read.is_duplicate or read.is_qcfail:
                    __read_qc_fail += 1
                    continue

                if int(read.mapping_quality) < self.quality_threshold:
                    __quality_trimmed += 1
                    continue

                # Check we have a CR:Z and UR:Z key:
                tags = dict(read.get_tags())
                if 'CB' in tags:
                    barcode = tags['CB']
                elif 'CR' in tags:
                    barcode = tags['CR']
                else: # No barcode
                    raise AssertionError('CB or CR tag not found!')
                    continue

                if whitelist and barcode not in whitelist:
                    # TODO: 1 bp mismatch recovery
                    __invalid_barcode_reads += 1
                    continue

                if UMIS:
                    if 'UB' in tags:
                        umi = f"{tags['UB']}-{barcode}" # UMI should be unique for both
                    elif 'UR' in tags:
                        umi = f"{tags['UR']}-{barcode}" # UMI should be unique for both
                    else:
                        raise AssertionError('UB or UR tag not found!')
                        continue
                else: #
                    umi = None # putting this here like this will ignore umis, and count all reads

                chrom = read.reference_name.replace('chr', '')
                if chrom not in self.genome.buckets: # Must be a valid chromosome
                    continue

                left = read.reference_start
                rite = read.reference_end
                loc_strand = '-' if read.is_reverse else '+'

                # Check we havne't seen this UMI/CB before:
                if not umi:
                    pass # Just skip;

                elif umi in umis: # umi/CB/chrom/strand was seen
                    if strand:
                        l = (chrom, loc_strand)
                    else:
                        l = (chrom, )

                    if UMIS and l in umis[umi]: # check we haven't seen this fragment;
                        __already_seen_umicb += 1
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
                    loc1_left = left
                    loc1_rite = left+1
                    loc2_left = rite-1
                    loc2_rite = rite # I just align the two edges, then I don't need to worry about split-reads, and I rely on the duplicate removal

                    # To get rid of the same gene twice
                    for buck in buckets_reqd:
                        if buck in self.genome.buckets[chrom]:
                            loc_ids.update(self.genome.buckets[chrom][buck]) # set = unique ids

                    for index in loc_ids:
                        locG_l = loc_lookups[index][0]
                        locG_r = loc_lookups[index][1]
                        # check strands

                        #print loc.qcollide(self.linearData[index]["loc"]), loc, self.linearData[index]["loc"]
                        #if loc1.qcollide(self.genome.linearData[index]["loc"]): # Any 1 bp overlap...
                        #    result.append(self.genome.linearData[index])

                        locG = self_genome_linearData[index]["loc"]

                        if loc1_rite >= locG_l and loc1_left <= locG_r:
                            result.append(self_genome_linearData[index])

                        if loc2_rite >= locG_l and loc2_left <= locG_r: # Any 1 bp overlap...
                            result.append(self_genome_linearData[index])

                    if result:
                        # We are going to add to something:
                        if barcode not in self.barcodes:
                            self.barcodes[barcode] = 0
                        self.barcodes[barcode] += 1
                        # do the annotation so that a read only gets counted to a TE if it does not hit a gene:
                        # This will currently allow 1 read to be counted twice if each edge is inside a different feature.
                        # Is that wrong, or a reasonable compromise?

                        types = set([i['type'] for i in result])
                        ensgs = set([(i['ensg'], i['strand']) for i in result]) # only count 1 read to 1 gene
                        if 'protein_coding' in types or 'lincRNA' in types or 'lncRNA' in types:
                            for e in ensgs:
                                # e[1] = strand Only collect gene results on the correct strand
                                if loc_strand != e[1]:
                                    continue

                                if ':' in ensgs: # A TE, skip it
                                    continue

                                if barcode not in final_results[e[0]]:
                                    final_results[e[0]][barcode] = set([])
                                final_results[e[0]][barcode].add(umi)

                        elif 'TE' in types:
                            for e in ensgs: # Not in any other RNA, so okay to count as a TE
                                if barcode not in final_results[e[0]]:
                                    final_results[e[0]][barcode] = set([])
                                final_results[e[0]][barcode].add(umi)
                        elif 'enhancer' in types:
                            for e in ensgs: # Not in any other RNA, so okay to count as a enhancer
                                if barcode not in final_results[e[0]]:
                                    final_results[e[0]][barcode] = set([])
                                final_results[e[0]][barcode].add(umi)
                        __read_assinged_to_gene += 1
                        #print()

        except StopIteration:
            pass # the last read

        sam.close()

        __total_rejected_reads = __already_seen_umicb + __quality_trimmed +__read_qc_fail
        __total_valid_reads = idx - __total_rejected_reads

        log.info('Processed {:,} SE reads'.format(idx))
        log.info('{:,} invalid barcode reads'.format(__invalid_barcode_reads))
        log.info('{:,} UMI-CB combinations were seen multiple times and removed'.format(__already_seen_umicb))
        log.info('{:,} Read quality is too low (<{})'.format(__quality_trimmed, self.quality_threshold))
        log.info('{:,} Reads QC failed'.format(__read_qc_fail))
        log.info('{:,} total reads rejected'.format(__total_rejected_reads))
        log.info('{:,} total valid reads'.format(__total_valid_reads))
        log.info('Assigned {:,} ({:.1f}%) valid reads to features'.format(__read_assinged_to_gene, ((__read_assinged_to_gene/__total_valid_reads) * 100.0))) # add per cents here;

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
            log.info(f'Keeping the best {maxcells:,} barcodes')
            barcodes_to_do = [i[0] for i in barcodes_to_do][0:maxcells]
        elif maxcells > len(self.barcodes):
            log.warning('Asked for {0:,} maxcells, but only {1:,} barcodes found'.format(maxcells, len(self.barcodes)))
            barcodes_to_do = [i[0] for i in barcodes_to_do]
        else:
            barcodes_to_do = [i[0] for i in barcodes_to_do]

        if '.tsv' not in out_filename:
            out_filename = f"{out_filename}.tsv"
        barcode_freq_filename = out_filename.replace('.tsv', '.barcode_freq.tsv')

        with open(barcode_freq_filename, 'w') as oh:
            for b in barcodes_to_do:
                oh.write('{0}\t{1}\n'.format(b, self.barcodes[b]))
        log.info('Saving barcode read frequency file to {0}'.format(barcode_freq_filename))

        with open(out_filename, 'w') as oh:
            oh.write('{0}\t{1}\n'.format('name', '\t'.join(result.keys())))

            for barcode in barcodes_to_do:
                counts = []
                for feature in result:
                    if barcode in result[feature]: # Stop defaultdict from densifying
                        counts.append(len(result[feature][barcode]))
                    else:
                        counts.append(0)
                oh.write('{0}\n'.format('\t'.join([barcode] + [str(c) for c in counts])))

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt\n")
        sys.exit(0)
