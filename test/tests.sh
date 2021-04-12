
# test bulk RNA-seq:
te_count -i data/bulk_rnaseq_PE.bam -o run_results/bulk_out_PE.tsv -g mm10 -m genes_tes
# Then md5 on the result:

# Test single cell path:
te_count -i data/single_cell_rnaseq.bam -g mm10 --se --sc -m genes_tes -o run_results/single_cell_out.tsv  --maxcells 3
te_count -i data/single_cell_rnaseq.bam -g mm10 --se --sc -m genes_tes -o run_results/single_cell_strand_out.tsv -w barcodes/version1.txt --maxcells 3 --strand

# enh bulk
te_count -i data/bulk_rnaseq_PE.bam -o run_results/bulk_out_PE-enh.tsv -g mm10 -m enhancers

# enh single cell path:
te_count -i data/single_cell_rnaseq.bam -g mm10 --sc -m enhancers -o run_results/single_cell_out-enh.tsv  --maxcells 3

#md5sum single_cell_out-enh-expected.tsv 
#md5sum md5sum single_cell_out-enh.tsv  
#match='1ae7dd557cb70bb83ccea84c84ec390c'


