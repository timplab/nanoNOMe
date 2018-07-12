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
  wget ftp://ftp.ensembl.org/pub/release-92/gtf/homo_sapiens/Homo_sapiens.GRCh38.92.gtf.gz \
    -O $db
# cgi
db=$dbdir/hg38_cgi.txt.gz
[ -e $db ]||\
  wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExtUnmasked.txt.gz \
    -O $db

# gm12878
cell=gm12878
samp=$(echo "$cell" | tr "[a-z]" "[A-Z]")
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
