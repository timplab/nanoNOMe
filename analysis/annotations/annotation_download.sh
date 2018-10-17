#!/bin/bash
# this is a script to download annotations
# you must provide root as an argument
root="$1"
echo "downloading annotations"

# hg38 annotations
dbdir=$root/hg38
[ -e $dbdir ]||mkdir $dbdir
# repeats
db=$dbdir/hg38_repeats.txt.gz
[ -e $db ]||\
  wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/nestedRepeats.txt.gz \
    -O $db
# genes
db=$dbdir/hg38_genes.gtf.gz
[ -e $db ]||\
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz \
    -O $db
# cgi
db=$dbdir/hg38_cgi.txt.gz
[ -e $db ]||\
  wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExtUnmasked.txt.gz \
    -O $db

# enhancers
db=$dbdir/hg19_enhancers.txt.gz
[ -e $db ]||\
  wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/vistaEnhancers.txt.gz \
    -O $db


# gm12878
cell=gm12878
samp=$(echo "$cell" | tr "[a-z]" "[A-Z]")
dbdir=$root/$cell
[ -e $dbdir ]||mkdir $dbdir
#celldir=$root/$cell
#dbdir=$celldir/$dbname
#[ -e $dbdir ]||mkdir -p $dbdir
#
#pre=$dbdir/${samp}_$2
#if [ "$2" == "mnase" ];then
#  #wget https://www.encodeproject.org/files/ENCFF000VME/@@download/ENCFF000VME.bigWig -O $dbdir/${samp}_mnase.bigWig
#  bigWigToBedGraph $pre.bigWig $pre.bedGraph
#fi
#
#if [ "$2" == "rnaseq" ];then
#  rep1=$pre.1.tsv
#  rep2=$pre.2.tsv
#  wget https://www.encodeproject.org/files/ENCFF212CQQ/@@download/ENCFF212CQQ.tsv -O $rep1
#  wget https://www.encodeproject.org/files/ENCFF350QZU/@@download/ENCFF350QZU.tsv -O $rep2
#fi
#
# enhancers
# taken from enhancer atlas
dir=$dbdir/enhancer
[ -e $dir ]||mkdir $dir
db=$dir/${samp}_enhancer_promoter_hg19.txt
[ -e $db ]||\
  wget http://enhanceratlas.org/data/AllEPs/human/GM12878_EP.txt \
  -O $db

