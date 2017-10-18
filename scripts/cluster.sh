#!/bin/bash
REF=$1
fastq1=$2
fastq2=$3
prefix=$4
dir=$5
hmms=$6
CPUS=8
cd $dir
if [ ! -f ${prefix}_${REF}.bam ]; then
  echo "[Running BWA]" `date "+%s"`
  echo "bwa mem -t $CPUS -M $WORK/data/references/hg19/${REF}.fa $fastq1 $fastq2 | samtools view -bS - > ${prefix}_${REF}.bam" 
  bwa mem -t $CPUS -M $WORK/data/references/hg19/${REF}.fa $fastq1 $fastq2 | samtools view -bS - > temp.bam
  mv temp.bam ${prefix}_${REF}.bam;
  echo "[Finished BWA]" `date "+%s"` 
fi

if [ ! -f ${prefix}_${REF}.unknown.bam ] ||  [ ! -s ${prefix}_${REF}.unknown.bam ]; then
  echo "[Find trans reads]" `date "+%s"` 
  python $WORK/bin/python/get_trans_new.py --unknown ${prefix}_${REF}.unknown.bam --data ${prefix}_${REF}.bam --trans ${prefix}_${REF}.trans.bam --viral ${prefix}_${REF}.viral.bam
fi
if [ ! -f test/temp/reduced.csv ]; then 
  echo "[Running HMMs]" `date "+%s"` 
  mkdir -p test
  $HOME/localperl/bin/perl $WORK/bin/perl/identify_viruses.pl -b ${prefix}_${REF}.unknown.bam -d test -H $hmms > test/hmms.log 2>&1
  echo "[Finished HMMs]" `date "+%s"`   
fi

python $WORK/bin/python/merge_viral_reads.py --unknown ${prefix}_${REF}.unknown.bam --trans ${prefix}_${REF}.trans.bam --reduced test/temp/reduced.csv --map test/temp/unmapped.map --output ${prefix}_${REF}.fixed.trans.bam

samtools sort -m 2G -@ $CPUS ${prefix}_${REF}.fixed.trans.bam > ${prefix}_${REF}.fixed.trans.cs.bam
samtools sort -m 2G -@ $CPUS ${prefix}_${REF}.viral.bam > ${prefix}_${REF}.viral.cs.bam
samtools index ${prefix}_${REF}.fixed.trans.cs.bam
samtools index ${prefix}_${REF}.viral.cs.bam
echo "[Finished trans reads]" `date "+%s"` 
echo "[Cluster reads]" `date "+%s"` 
#python $WORK/bin/python/cluster_trans_new.py --data ${prefix}_${REF}.fixed.trans.cs.bam --output ${prefix}_${REF}.clusters.txt
python /pedigree2/projects/namphuon/data/cancer_viral_integration/src/cluster_trans_new.py --data ${prefix}_${REF}.fixed.trans.cs.bam --output ${prefix}_${REF}.clusters.txt
echo "[Finished Cluster reads]" `date "+%s"` 

