#!/bin/bash
env=aws
root=/shared/data/nomeseq
outroot=$root
batch=/shared/Code/ilee/slurm/batchcommand.scr
bedscr=/shared/Code/ilee/nanopolish/script/mtsv2bedGraph.py

for samp in mcf10a mcf7 mdamb231;do
  outdir=$root/$samp/ngmlr
  mbeddir=$outdir/mbed
  logdir=$outdir/log
  for mod in cpg gpc;do
    mcalldir=$outdir/mcall-$mod
    for rep in 1 2 3;do
      dir=$(find $mcalldir/* -maxdepth 0 -name "*rep$rep*")
      base=$(basename "$dir")
      bed=$mbeddir/$base.$mod.meth.bed.gz
      tsvs=$(find $dir -name "*tsv" -type f)
      subcom="cat $tsvs |\
        python $bedscr -m $mod |\
        sort -T $mbeddir -k1,1 -k2,2n |\
        bgzip > $bed && tabix -p bed $bed"
      log=$logdir/$base.mbed.log
      com="sbatch -e $log -o $log -J mbed $batch $subcom"
      echo $com
      $com

    done
  done
done


