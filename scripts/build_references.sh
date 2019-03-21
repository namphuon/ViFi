#!/bin/bash
INPUT_FASTA=$1
OUTPUT_DIR=$2
PREFIX=$3

#Build alignment/tree
run_pasta.py --max-mem-mb=4000 -d dna -j $PREFIX -i $INPUT_FASTA  -o $OUTPUT_DIR

OUT_ALN=`grep "Writing resulting alignment" $PREFIX.out.txt | awk '{print $7}'`
OUT_TREE=`grep "Writing resulting tree" $PREFIX.out.txt | awk '{print $7}'`

#Decompose alignment/tree into HMMs
python $VIFI_DIR/scripts/build_hmms.py --tree_file $OUT_TREE --alignment_file $OUT_ALN --prefix $PREFIX --output_dir $OUTPUT_DIR

#Build HMM list
ls $OUTPUT_DIR/*.hmmbuild > $OUTPUT_DIR/hmm_list.txt
