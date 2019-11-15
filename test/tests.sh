
# test bulk RNA-seq:
te_count -i bulk_rnaseq_PE.bam -o bulk_out_PE.tsv -g mm10 -m genes_tes
# Then md5 on th result:

# Test single cell path:
te_count -i single_cell_rnaseq.bam -g mm10 --sc -m genes_tes -o single_cell_out.tsv  --maxcells 3



