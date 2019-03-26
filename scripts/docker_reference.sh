#!/bin/bash
INPUT_DIR=$1
READ1=$2
READ2=$3
OUTPUT_DIR=$4
CPUS=${5:-1}

docker run -e CPUS=$CPUS -v $REFERENCE_REPO:/home/repo/data -v $INPUT_DIR:/home/fastq/ -e READ1=$READ1 -e READ2=$READ2 -v $AA_DATA_REPO:/home/data_repo/ -v $OUTPUT_DIR:/home/output/ vifi:latest sh /home/vifi.sh

