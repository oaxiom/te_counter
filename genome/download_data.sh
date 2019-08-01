
# May not be correct or latest versions, but checked on 21/1/2019.

# gencode
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M20/gencode.vM20.annotation.gtf.gz

wget -c -O hg38_rmsk.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
wget -c -O mm10_rmsk.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz

wget ftp://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/chromInfo.txt.gz -O mm10.chromSizes.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/chromInfo.txt.gz -O hg38.chromSizes.gz

gunzip -c mm10.chromSizes.gz | grep -v -E 'random|chrUn|chrM'  >mm10.chromSizes.clean
gunzip -c hg38.chromSizes.gz | grep -v -E 'random|chrUn|chrM|_alt|_fix'  >hg38.chromSizes.clean