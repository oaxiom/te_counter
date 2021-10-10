# te_counter
A simple tool to count genes/TEs for reads from a BAM file

# INSTALL

No easy install, it's pure python, so just stick it on your PATH.

Requires pysam.

# TUTORIAL

Build the indeces like this:

te_genome -g hg38 -m gene_tes

Supported genomes are mm10 and hg38

You can build custom gtfs with --gtf (see help below)

Count BAM files like this:

Go to te_counter/test and you can execute the test scripts for various configurations:

Paired-end bulk RNA-seq:
te_count -i data/bulk_rnaseq_PE.bam -o run_results/bulk_out_PE.tsv -g mm10 -m genes_tes

# Single-end bulk RNA-seq (uses the paired-end BAM, pretend it's single-end)
te_count --se -i data/bulk_rnaseq_PE.bam -o run_results/bulk_out_SE.tsv -g mm10 -m genes_tes

# 10x-style single cell data, expects CB, UB SAM flags
te_count -i data/single_cell_rnaseq.bam -g mm10 --se --sc -m genes_tes -o run_results/single_cell_out.tsv  --maxcells 3
te_count -i data/single_cell_rnaseq.bam -g mm10 --se --sc -m genes_tes -o run_results/single_cell_strand_out.tsv -w barcodes/version1.txt --maxcells 3 --strand

# HELP

te_genome -h
usage: te_genome [-h] [--gtf GTF] [-f] -m MODE -g GENOME

Builds the indices for te_count

required arguments:
  -m MODE, --mode MODE  Type of annotation to use, valid modes: {'enhancers', 'genes_tes', 'custom'}
  -g GENOME, --genome GENOME
                        Genome assembly to use, valid genomes: {'mm10', 'hg38'}

Example usage: te_genome -g genome -m mode

te_count -h

usage: te_count [-h] [--se] [--sc] [--noumi] [--strand] [-q QUAL] [--maxcells MAXCELLS] [-w W] -i INBAM -o OUTTSV -g GENOME -m MODE

Counts up the number of reads that overlap some set of gene/TE or other features

required arguments:
  -i INBAM, --inbam INBAM
                        the BAM alignment file containing the reads
  -o OUTTSV, --outtsv OUTTSV
                        the TSV file to save the genes and count data to
  -g GENOME, --genome GENOME
                        Genome assembly to use, valid genomes: {'mm10', 'hg38'}
  -m MODE, --mode MODE  Type of annotation to use, valid modes: {'custom', 'enhancers', 'genes_tes'}

optional arguments:
  -h, --help            show this help message and exit
  --se                  Set mode to SE (single-end) mode, default is paired-end mode
  --sc                  Set te_count to single-cell mode, default is bulk. Expects an CR:Z flag in the BAM/SAM file. Make sure you set --se
                        appropriately
  --noumi               If --sc is set, but this single cell data has no UMI in a UR:Z tag, set this switch
  --strand              Use the strand information for the reads (protocol is strand-specific), default=False
  -q QUAL, --qual QUAL  q threshold for qulait filtering, default=20
  --maxcells MAXCELLS   keep at most maxcells with the most reads, default=10,000
  -w W                  A whitelist of barcodes. Becomes a requried argument of --sc

Example usage: te_count -i in.bam -o out.bam -g genome -m genes_tes

