#!/bin/bash
root=/shared/Data
tardir=$root/tar
srcroot=/shared/Code
ontpath=$srcroot/ilee/oxford
ontslurm=$ontpath/slurm
logroot=$root/log
for samp in `find $tardir -maxdepth 1 -name "*_2.*tgz"`
do
  base=$(basename "$samp")
  base=${base%%.*}
  echo $base
  methdir=$root/$base/methcall
  for calldir in `find $methdir -maxdepth 1 -name "*-*" -type d`;do
    b=$(basename "$calldir")
    mod=${b%-*}
    logdir=$logroot/methfreq/$base/$b
    echo $b
    [ -e $logdir ]||mkdir -p $logdir
    if [ 1 -eq 1 ];then
      for tsv in `find $calldir -maxdepth 1 -name "*meth.tsv"`;do
        x=$(basename "$tsv")
        log=$logdir/$x.methfreq.log
        sbatch -e $log -o $log \
          $srcroot/ilee/dnamods/slurm/methfreq.scr -i $tsv -s $srcroot -m $mod
      done
    fi
    if [ 1 -eq 1 ];then
      freqlist=`find $calldir -maxdepth 1 -name "*freq.tsv"`
      log=$logdir/$base.$b.merge.log
      sbatch -e $log -o $log \
        $srcroot/ilee/dnamods/slurm/merge_methfreq.scr -i "$freqlist" -s $srcroot -o $base.$b.freq.pooled.tsv
    fi
  done
done


