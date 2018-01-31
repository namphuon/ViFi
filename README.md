# ViFi

ViFi is a tool for detecting viral integration and fusion mRNA sequences from Next Generation Sequencing data.  Unlike standard approaches that use reference-based read mapping for identification of viral reads, ViFi uses both reference-based read mapping and a phylogenetic-based approach to identify viral reads.  ViFi also incorporates mappability scores of the reads to filter out false positive integration detection.  The end result is a tool that can accurately and precisely detect integrated viruses, even if the viruses are highly mutated or novel strains.

ViFi is currently in alpha testing, is is constantly undergoing revisions.  High on the priority list is an easier installation process, as well as improve user interface.  Please report any problems/bugs to Nam Nguyen (ndn006@eng.ucsd.edu) so that ViFi can be improved and problems can be quickly corrected.  

#=============================================================================================================================================================================================

# Installation:

# ViFi download (if you have not already cloned this source code):
git clone https://github.com/namphuon/ViFi.git

# Dependencies:
## 1) Python 2.7
## 2) Pysam verion 0.9.0 or higher (https://github.com/pysam-developers/pysam):
sudo pip install pysam
## 3) Samtools 1.3.1 or higher (www.htslib.org/)
sudo apt-get install samtools
## 4) BWA 0.7.15 or higher (bio-bwa.sourceforge.net/)
sudo apt-get install bwa
## 5) Install HMMER v3.1b2 and have it on the path (http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-macosx-intel.tar.gz)


# Set the ViFi directory
echo export VIFI_DIR=/path/to/ViFi >> ~/.bashrc

# Data repositories:
## Download the data repositories. While we include some annotations, we are unable to host some large files in the git repository.
## These may be downloaded from https://drive.google.com/open?id=0ByYcg0axX7udUDRxcTdZZkg0X1k. Thanks to Peter Ulz for noticing incorrect link earlier.
tar zxf data_repo.tar.gz
echo export AA_DATA_REPO=$PWD/data_repo >> ~/.bashrc
source ~/.bashrc

## Download the HMM models from https://drive.google.com/open?id=0Bzp6XgpBhhghSTNMd3RWS2VsVXM 
unzip data.zip
echo export REFERENCE_REPO=$PWD/data/ >> ~/.bashrc

## For viral family of interest, create BWA index.  Example for HPV is given below
cat $AA_DATA_REPO//hg19/hg19full.fa $REFERENCE_REPO/hpv/hpv.unaligned.fas > $REFERENCE_REPO/hpv/hg19_hpv.fas
bwa index $REFERENCE_REPO/hpv/hg19_hpv.fas

## Running ViFi
python run_vifi.py -f <input_R1.fq.gz> -r <input_R2.fq.gz> -o <output_dir>

The ViFi manuscript is currently under review.

## Dockerized ViFi

We have also created a dockerized version of ViFi to enable easier time running.  The docker version of ViFi can be obtained
by installing Docker (https://www.docker.com/), and running the following command:
docker pull us.gcr.io/aa-test-175718/vifi

To run the dockerized version of ViFi, run the following script in the ViFi scripts directory:

docker_vifi.sh <READ1> <READ2> <OUTPUT>

