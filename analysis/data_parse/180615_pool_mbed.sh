#!/bin/bash
env=aws
root=/shared/data/nomeseq
outroot=$root
pooldir=$root/pooled/methylation
batch=/shared/Code/ilee/slurm/batchcommand.scr

for samp in mcf10a mcf7 mdamb231;do
  outdir=$root/$samp/ngmlr
  mbeddir=$outdir/mbed
  logdir=$outdir/log
  lab=$(echo $samp | awk '{ print toupper($1) }')
  echo $lab
  for mod in gpc;do
    tmpdir=$pooldir/tmp/$samp/$mod
    [ -e $tmpdir ]||mkdir -p $tmpdir
    out=$pooldir/$lab.$mod.pooled.meth.bed.gz
    beds=$(find $mbeddir -name "*$mod.meth.bed.gz")
    subcom="zcat $beds |\
        sort -T $tmpdir -k1,1 -k2,2n |\
        bgzip > $out && tabix -p bed $out"
    log=$logdir/$samp.bedmerge.log
    com="sbatch -e $log -o $log -J mbedmerge $batch $subcom"
    echo $com
    $com

  done
done


