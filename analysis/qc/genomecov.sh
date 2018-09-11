#d!/bin/bash
root=/dilithium/Data/Nanopore/projects/nomeseq/analysis/pooled/bed
beds=$(find $root -name "*bed" -type f)
qc="../../util/bedQC.py"

if [ "$1" == "qc" ];then
  for bed in $beds;do
    base=$(basename $bed)
    ext=${base#*.}
    samp=${base%%.*}
    samps="${samps} $samp"
  done
  parallel "python $qc $root/{}.$ext > $root/{}.qc.txt" ::: $samps
fi

if [ "$1" == "summarize" ];then
  cells="GM12878 MCF10A MCF7 MDAMB231"
  sumout=$root/bedqc.txt
  out=$(echo "readnum totalbp n50 q1 mean q3 max" | tr " " "\t")
  for cell in $cells;do
    qcnum=$(cat $root/$cell.qc.txt)
    sum=$(echo $cell $qcnum | tr " " "\t")
    out="$out;${sum}"
  done
  echo $out | tr ";" "\n" > $sumout
fi

