#!/bin/bash
root=/kyber/Data/Nanopore/projects/nanonome/analysis/data
repeats=$root/hg38/hg38_repeats.bed
regs=$root/hg38/hg38_repeats.multiple.bed

if [ ! -e $regs ]; then
  reps=$(cut -f4 $repeats | sort | uniq -c | sort -n | awk '$1>1000{ print $2 }' | tr "\n" "|")
  reps=${reps%"|"}
  grep -E "$reps" $repeats > $regs
fi
