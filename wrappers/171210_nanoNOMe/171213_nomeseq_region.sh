#!/bin/bash

dir=/shared/Data/171210_gpc_nome_mda231_200U*/methcall
#reg=igf2
#chr=chr11
#st=2120500
#end=2139500
#reg=gapdhtss
#chr=chr12
#st=6531000
#end=6538000
reg=wnt7bTSS
chr=chr22
st=45974000
end=45980000


for mod in cpg gpc;do
  echo $mod
  flist=`find $dir/${mod}-hg38 -name "*meth.tsv"`
  cat $flist | \
    awk -v chr=$chr -v st=$st -v end=$end \
    '{ if($1==chr&&$2>st&&$3<end) print }' > $dir/dmr.$reg.$mod.tsv
done

#cat *meth.tsv | awk '{ if($1=="chr9"&&$2>63815000&&$3<63819000) print }' > dmr.gpc.tsv
