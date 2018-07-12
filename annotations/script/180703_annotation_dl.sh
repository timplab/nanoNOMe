#!/bin/bash
root="/home/isac/Dropbox/Data/nome-seq/db"
cell=$1
dbname=$2
samp=$(echo "$cell" | tr "[a-z]" "[A-Z]")
echo $samp
celldir=$root/$cell
dbdir=$celldir/$dbname
[ -e $dbdir ]||mkdir -p $dbdir

pre=$dbdir/${samp}_$2
if [ "$2" == "mnase" ];then
  #wget https://www.encodeproject.org/files/ENCFF000VME/@@download/ENCFF000VME.bigWig -O $dbdir/${samp}_mnase.bigWig
  bigWigToBedGraph $pre.bigWig $pre.bedGraph
fi

if [ "$2" == "rnaseq" ];then
  rep1=$pre.1.tsv
  rep2=$pre.2.tsv
  wget https://www.encodeproject.org/files/ENCFF212CQQ/@@download/ENCFF212CQQ.tsv -O $rep1
  wget https://www.encodeproject.org/files/ENCFF350QZU/@@download/ENCFF350QZU.tsv -O $rep2
fi

