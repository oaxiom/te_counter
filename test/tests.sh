
mkdir run_results

# test bulk RNA-seq:
te_count -i data/bulk_rnaseq_PE.bam -o run_results/bulk_out_PE.tsv -g mm10 -m genes_tes
# Then md5 on the result:

# Yes, the BAM is PE, but you can just silently treat as SE
te_count --se -i data/bulk_rnaseq_PE.bam -o run_results/bulk_out_SE.tsv -g mm10 -m genes_tes

# Test single cell path:
te_count -i data/single_cell_rnaseq.bam -w barcodes/version1.txt -g mm10 --se --sc -m genes_tes -o run_results/single_cell_out.tsv  --maxcells 3
te_count -i data/single_cell_rnaseq.bam -w barcodes/version1.txt -g mm10 --se --sc -m genes_tes -o run_results/single_cell_strand_out.tsv -w barcodes/version1.txt --maxcells 3 --strand

# enh bulk
te_count -i data/bulk_rnaseq_PE.bam -o run_results/bulk_out_PE-enh.tsv -g mm10 -m enhancers

# enh single cell path:
te_count -i data/single_cell_rnaseq.bam -w barcodes/version1.txt -g mm10 --sc -m enhancers -o run_results/single_cell_out-enh.tsv  --maxcells 3

#md5sum single_cell_out-enh-expected.tsv 
#md5sum md5sum single_cell_out-enh.tsv  
#match='1ae7dd557cb70bb83ccea84c84ec390c'

# snrnps
te_count -i data/bulk_rnaseq_PE.bam -o run_results/bulk_snrnpsiPE.tsv -g hg38 -m snrnps
te_count --se -i data/bulk_rnaseq_PE.bam -o run_results/bulk_snrnpsSE.tsv -g hg38 -m snrnps

# Scalene optimization path:
python3 -m scalene main_for_profile.py
