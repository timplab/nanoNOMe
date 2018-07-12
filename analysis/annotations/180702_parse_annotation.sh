#!/bin/bash
root="/dilithium/Data/Nanopore/projects/nomeseq/database"
while :
do
  case "$1" in 
    -c | --cell )
      cell=$( echo "$2" | tr "[A-Z]" "[a-z]" )
      shift 2
      ;;
    -d | --database )
      dbname=$2
      shift 2
      ;;
    * ) break
      ;;
  esac
done
if [ -z "$cell" ];then
  echo "cell name must be given by -c"
  exit
fi
if [ -z "$dbname" ];then
  echo "database msut be given by -d"
  exit
fi
task=$1
lab=$(echo "$cell" | tr "[a-z]" "[A-Z]" )
dbdir=$root/$cell/$dbname
[ -e $dbdir ]||mkdir -p $dbdir
pre=$dbdir/${lab}_$dbname

parser="../../util/bed_parser.py"

if [ "$dbname" == "dnase" ];then
  db=$dbdir/"ENCFF598KWZ.bed.gz"
elif [ "$dbname" == "atac" ];then
  db=$dbdir/"GSE47753_GM12878_ATACseq_50k_AllReps_ZINBA_pp08.hg38.bed"
elif [ "$dbname" == "ctcf" ];then
  db=$dbdir/"GM12878_CTCF.bed"
fi

echo $db
center=$pre.center.bed
if [ "$task" == "getcenter" ];then
  $parser center -b $db -o $center
fi
size=2000
side=$(($size/2))
region=$pre.${size}bp.bed
if [ "$task" == "getregion" ];then
  $parser region -b $center -u $side -d $side |\
    sort -k1,1 -k2,2n > $region
fi
