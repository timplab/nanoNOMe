#!/bin/bash
root=/mithril/Data/NGS/Reference/human_annotations
name=cpgIslandExtUnmasked
pre=$root/$name
dl=$pre.txt.gz
bed=$pre.bed
tssname=hg38.91.TSS.2kb
tssbed=$root/$tssname.bed
cgitss=$root/$tssname.CGI.bed

# first download the database
[ -e $dl ]||\
  wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExtUnmasked.txt.gz \
  -O $dl

# unzip and convert to bed
[ -e $bed ]||\
  gunzip -c $dl |\
    awk 'OFS="\t"{ print $2,$3,$4,".",".","." }' \
    > $bed
# get cgi TSS
# total number of known genes : 58300
# first start off with complete overlap only => <1k regions
# half overlap : 6666
# any overlap : 21992 > let's go with any overlap
if [ "$1" == "cgitss" ];then
  bedtools intersect -a $tssbed -b $bed |\
    sort -k1,1 -k2,2n > $cgitss
fi
