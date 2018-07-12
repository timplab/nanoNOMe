#!/bin/bash
root=/shared/Data/analysis
cell=gm12878
dir=$root/$cell/ngmlr
beddir=$dir/mbed
metadata=$root/nome-seq_data.csv
scriptdir=/shared/Code/ilee/nanopolish/script
freqdir=$dir/meth_frequency
[ -e $freqdir ]||mkdir $freqdir
logdir=$dir/log/pool
[ -e $logdir ]|| mkdir $logdir
batchscr=/shared/Code/ilee/slurm/batchcommand.scr
pooldir=$dir/pooled

for mod in cpg gpc;do
  beds=`readlink -f $beddir/*$mod.meth.bedGraph`
  poolbed=$pooldir/GM12878.$mod.meth.bedGraph
  mergecom="cat $beds | sort -T $pooldir -k1,1 -k2,2n > $poolbed"
  log=$logdir/GM12878.$mod.meth.pool.log
  com="sbatch -J merge -e $log -o $log $batchscr $mergecom"
  echo $com
  $com
done

