#!/bin/bash
# parse repeats
# path to the UCSC nestedrepeats.txt.gz provided as an arugment
echo "parsing repeats"
db="$1"
prefix=${db%%.*}

# convert to bed
bed=$prefix.bed
[ -e $bed ]||\
  gunzip -c $db |\
  awk 'OFS="\t"{ print $2,$3,$4,$5,$6,$7,$16 }' |\
  sort -k1,1 -k2,2n > $bed

# subset by repeat type
repeat_types="LINE SINE DNA LTR Retroposon"
echo "types of repeats being parsed : $repeat_types"
for repeat in $repeat_types ;do
  subbed=${prefix}_$repeat.bed
  [ -e $subbed ]||\
    grep -v ? $bed | grep $repeat > $subbed
done
# special case : LTR and Retroposons can be merged
subbed=${prefix}_RTP.bed
[ -e $subbed ]||\
  cat ${prefix}_LTR.bed ${prefix}_Retroposon.bed |\
  sort -k1,1 -k2,2n > $subbed

# LINE : subset larger regions
#awk 'OFS="\t"{ if($3-$2>=5000) print }' ${prefix}_LINE.bed > ${prefix}_LINE_long.bed
