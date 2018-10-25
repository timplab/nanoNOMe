#!/bin/bash 
set -eo pipefail

ROOT=/dilithium/Data/Nanopore/projects/nomeseq/analysis
STEP="all"
while :
do
  case "$1" in
    -p | --path ) # root path
      ROOT="$2"
      shift 2
      ;;
    -s | --step ) # step to perform
      STEP="$2"
      shift 2
      ;;
    * ) break
      ;;
  esac
done

echo "analysis root directory : $ROOT"
[ -e $ROOT/tmp ]||mkdir $ROOT/tmp
echo "performing $STEP step(s)"
CODEROOT=$(readlink -f $PWD/analysis)

############################################################
#
# Step 0. Download data and install software
# 
############################################################

if [[ $STEP =~ all|step0|hg38Annotation ]];then
  DBDIR=$ROOT/annotations/hg38
  [ -e $DBDIR ]||mkdir -p $DBDIR
  # genes
  genedb=$DBDIR/hg38_genes.gtf.gz
	[ -e $genedb ]||\
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz \
      -O $genedb
  # cgi
  cgidb=$DBDIR/hg38_cgi.txt.gz
	[ -e $cgidb ]||\
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExtUnmasked.txt.gz \
      -O $cgidb
  # repetitive elements
	repdb=$DBDIR/hg38_repeats.txt.gz
	[ -e $repdb ]||\
		wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/nestedRepeats.txt.gz \
			-O $repdb
  # liftover
  chain=$DBDIR/hg19ToHg38.over.chain
	if [ ! -e $chain ];then
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz \
      -O $chain.gz
    gunzip $chain.gz
  fi
  # parsing -
  # genes
  $CODEROOT/annotations/parse_genes.sh $genedb
  # repeats
  $CODEROOT/annotations/parse_repeats.sh $repdb
  # cgi
  gunzip -c $cgidb | awk 'OFS="\t"{ print $2,$3,$4 }' > $DBDIR/hg38_cgi.bed
  # shuffled bed for random sampling
  cp annotations/hg38_shuffle.bed $DBDIR
  
fi

if [[ $STEP =~ all|step0|gmAnnotation ]];then
  DBDIR=$ROOT/annotations/gm12878
  [ -e $DBDIR ]||mkdir -p $DBDIR
  # ctcf
  # download chip data
  CHIP=$DBDIR/GM12878_CTCF_ChIP.bed
  wget https://www.encodeproject.org/files/ENCFF356LIU/@@download/ENCFF356LIU.bed.gz \
    -O $CHIP.gz && gunzip $CHIP.gz
  awk 'OFS="\t"{ print $1,$4-1,$5,".",".",$7 }' annotations/CTCF_bindingsites.txt \
    > annotations/CTCF_bindingsites_hg19.bed
  CTCF=annotations/CTCF_bindingsites_hg38.bed
  # first liftover
  liftOver annotations/CTCF_bindingsites_hg19.bed \
    $ROOT/annotations/hg38/hg19ToHg38.over.chain \
    $CTCF $CTCF.unmapped
  # parse to get GM12878
  $CODEROOT/annotations/parse_GM12878_CTCF.sh $DBDIR $CTCF

  #DNAse
  DNASE=$DBDIR/GM12878_DNAse.bed
  wget https://www.encodeproject.org/files/ENCFF590NWM/@@download/ENCFF590NWM.bed.gz \
    -O $DNASE.gz && gunzip $DNASE.gz
  $CODEROOT/annotations/parse_GM12878_DNAse.sh $DBDIR $DNASE

fi

if [[ $STEP =~ all|step0|gmData ]];then
  DBDIR=$ROOT/data/gm12878
  [ -e $DBDIR ]||mkdir -p $DBDIR
  # mnase
  PRE=$DBDIR/GM12878_MNase
  [ -e $PRE.bigWig ]||\
    wget https://www.encodeproject.org/files/ENCFF000VME/@@download/ENCFF000VME.bigWig \
    -O $PRE.bigWig
  if [ ! -e $PRE.bedGraph.gz ];then
    bigWigToBedGraph $PRE.bigWig $PRE.bed 
    bgzip $PRE.bedGraph 
    tabix -p bed $PRE.bedGraph.gz
  fi
  MNASE=${PRE}_hg38.bedGraph
  if [ ! -e $MNASE.gz ];then
    liftOver $PRE.bed.gz $ROOT/annotations/hg38/hg19ToHg38.over.chain \
      $MNASE $MNASE.unmapped
    sort -k1,1 -k2,2n -T $ROOT/tmp $MNASE.sorted && bgzip $MNASE.sorted
    mv $MNASE.sorted $MNASE.gz
    tabix -p bed $MNASE.gz
  fi
fi

if [[ $STEP =~ all|step0|bcanAnnotation ]];then
  DBDIR=$ROOT/annotations/breastcancer
  [ -e $DBDIR ]||mkdir -p $DBDIR
  hg38genes=$ROOT/annotations/hg38/hg38_genes.bed
  bcan=annotations/bcan_genes_CancerGeneCensus.txt
  names=$(awk 'NR>1{ print $1 }' $bcan)
  BCANBED=$DBDIR/bcan_genes.bed
  BCANPROM=$DBDIR/bcan_genes.promoter.bed
  grep -F "$names" $hg38genes > $BCANBED
  parser=$CODEROOT/../util/bed_parser.py 
  python $parser getstart -b $BCANBED |\
    python $parser region -u 2000 -d 3000  |\
    sort -k1,1 -k2,2n > $BCANPROM
  # repetitive regions near bcan genes
  rep=$ROOT/annotations/hg38/hg38_repeats.bed
  proxrep=$DBDIR/repeats_bcan_proximity.bed
  bedtools closest -D b -a $rep -b $BCANBED |\
    awk '{ if(sqrt($NF*$NF) <= 5000) print }' > $proxrep
  bodyrep=$DBDIR/repeats_bcan_body.bed
  bedtools intersect -wo -a $rep -b $BCANBED > $bodyrep
  cp annotations/bcan_10a_vs_231_repeats.bed $ROOT/annotations/breastcancer/
fi



############################################################
#
# Step 1. Basecalling, alignment, methylation calling,
#         and converting to methylation bed format
#
############################################################



############################################################
#
# Step 2. Figure 2 metaplots, expression plots, repeats
#
############################################################

# metaplots
if [[ $STEP =~ all|step2|metaplots ]];then
  echo "metaplots (ctcf,dnase-seq)"
  DBDIR=$ROOT/annotations/gm12878
  OUTDIR=$ROOT/plots/metaplots
  [ -e $OUTDIR ]||mkdir -p $OUTDIR
  # methylation
  CELLS="GM12878 GM12878_BSseq_1"
  DIR=$ROOT/pooled/methylation/mfreq_all
  for CELL in $CELLS ; do
    cpg=$DIR/$CELL.cpg.methfreq.txt.gz
    gpc=$DIR/$CELL.gpc.methfreq.txt.gz
    # CTCF
    Rscript script/nanonome_plots.R metaplotByDistance -c $cpg -g $gpc \
      -r $DBDIR/GM12878_CTCF_hg38.center.2000bp.bed \
      -o $OUTDIR/$CELL.CTCF.metaplot.pdf 2> /dev/null
    # DNAse
    Rscript script/nanonome_plots.R metaplotByDistance -c $cpg -g $gpc \
      -r $DBDIR/GM12878_DNAse_top.center.2000bp.bed \
      -o $OUTDIR/$CELL.DNAse.metaplot.pdf 2> /dev/null
  done
  # mnase
  analysis/frequency/180921_mnase_metaplot.sh $ROOT
fi

# expression plots
if [[ $STEP =~ all|step2|gm12878exp ]];then
  echo "GM12878 tss chromatin by expression plots"
  Rscript $CODEROOT/frequency/181015_gm12878_tss_chromatin_by_exp.R $ROOT 2> /dev/null
fi
# repeats
if [[ $STEP =~ all|step2|gm12878coverage ]];then
  echo "GM12878 BS-seq vs nanoNOMe coverage comparison"
  # first global plots
  Rscript $CODEROOT/frequency/180707_bsseqVSnano_coveragecomp.R $ROOT 2> /dev/null
fi
# specific regions, read-level
if [[ $STEP =~ all|step2|gm12878readlevel ]];then
  echo "GM12878 BS-seq vs nanoNOMe read-level comparison"
  $CODEROOT/igv/gm12878_igv_plotting.sh $ROOT
fi

############################################################
#
# Step 3. Figure 3 Breast cancer comparison - 
#         promoter region comparison,
#         enhancer promter chromatin state correlation (heatmap?)
#
############################################################

# bcan promoters, read-level
if [[ $STEP =~ all|step3|bcangenes ]];then
  echo "bresat cancer sample gene comparisons"
  $CODEROOT/igv/bcan_igv_plotting.sh $ROOT bcanpromoter
fi


############################################################
#
# Step 4. Figure 4 SV analysis - 
#         MDA-MB-231-exclusive SV falling in bcan gene with 
#         interesting chromatin pattern
#
############################################################

if [[ $STEP =~ all|step4|sv ]];then
  echo "bresat cancer sample sv comparisons"
  $CODEROOT/sv/181020_bcanSV.sh
#  $CODEROOT/frequeyncy/bcan_sv_mfreq_comparison.R
fi



