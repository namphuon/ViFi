#!/bin/bash
INPUT_FASTA=$1
PREFIX=$2
HMMS=$3

OUTPUT_DIR=$REFERENCE_REPO/$PREFIX
mkdir -p $OUTPUT_DIR
OUTPUT_DIR=`realpath $OUTPUT_DIR`

#Create hg19+virus fasta file
cat $AA_DATA_REPO//hg19/hg19full.fa $INPUT_FASTA > $REFERENCE_REPO/$PREFIX/hg19_${PREFIX}.fas

#Build index
docker run -v $REFERENCE_REPO/$PREFIX/:/home/$PREFIX/ docker.io/namphuon/vifi bwa index /home/$PREFIX/hg19_${PREFIX}.fas

