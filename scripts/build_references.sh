#!/bin/bash
INPUT_FASTA=$1
OUTPUT_DIR=$2
PREFIX=$3

mkdir -p $OUTPUT_DIR
OUTPUT_DIR=`realpath $OUTPUT_DIR`
INPUT_DIR=`dirname $INPUT_FASTA`
if [ "$INPUT_DIR" == "." ]
then
  INPUT_DIR=`pwd`
fi 
INPUT_NAME=`basename $INPUT_FASTA`

#Pull latest version of PASTA
docker pull smirarab/pasta

#Build alignment/tree, use 4GB for alignment, via Docker
docker run -v $INPUT_DIR/:/data/ -v $OUTPUT_DIR/:/output/ smirarab/pasta run_pasta.py  --max-mem-mb=4000 -d dna -j $PREFIX -i $INPUT_NAME  -o /output/

cd $OUTPUT_DIR

OUT_ALN=`grep "Writing resulting alignment" $PREFIX.out.txt | awk '{print $7}'`
OUT_ALN=`basename $OUT_ALN`
OUT_TREE=`grep "Writing resulting tree" $PREFIX.out.txt | awk '{print $7}'`
OUT_TREE=`basename $OUT_TREE`

#Decompose alignment/tree into HMMs via Docker
docker run -v $OUTPUT_DIR/:/output/ docker.io/namphuon/vifi python "scripts/build_hmms.py" --tree_file /output/$OUT_TREE --alignment_file /output/$OUT_ALN --prefix $PREFIX --output_dir /output/

#Build HMM list
ls $OUTPUT_DIR/*.hmmbuild > $OUTPUT_DIR/hmm_list.txt
