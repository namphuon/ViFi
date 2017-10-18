# ViFi

ViFi is a tool for detecting viral integration and fusion mRNA sequences from Next Generation Sequencing data.  Unlike standard approaches that use reference-based read mapping for identification of viral reads, ViFi uses both reference-based read mapping and a phylogenetic-based approach to identify viral reads.  ViFi also incorporates mappability scores of the reads to filter out false positive integration detection.  The end result is a tool that can accurately and precisely detect integrated viruses, even if the viruses are highly mutated or novel strains.

ViFi is currently in alpha testing, so please report any problems/bugs to Nam Nguyen (ndn006@eng.ucsd.edu) so that it can be quickly corrected.  

#=============================================================================================================================================================================================

# Installation:

# ViFi download (if you have not already cloned this source code):
git clone https://github.com/namphuon/ViFi.git

# Dependencies:
## 1) Python 2.7
## 2) Pysam verion 0.9.0 or higher (https://github.com/pysam-developers/pysam):
sudo pip install pysam


# Data repositories:
## Download the data repositories. While we include some annotations, we are unable to host some large files in the git repository.
## These may be downloaded from https://drive.google.com/open?id=0ByYcg0axX7udUDRxcTdZZkg0X1k. Thanks to Peter Ulz for noticing incorrect link earlier.
tar zxf data_repo.tar.gz
echo export AA_DATA_REPO=$PWD/data_repo >> ~/.bashrc
source ~/.bashrc

## Download the HMM models from https://drive.google.com/file/d/0Bzp6XgpBhhghRWRsczVOZy1yQ2c/view?usp=sharing
unzip data.zip
echo export HMM_REPO=$PWD/data/ >> ~/.bashrc

The ViFi manuscript is currently under review.


docker run --rm -it ubuntu
apt-get update
apt-get -y install build-essential checkinstall
apt-get -y install libreadline-gplv2-dev libncursesw5-dev libssl-dev libsqlite3-dev tk-dev libgdbm-dev libc6-dev libbz2-dev
apt-get -y install python2.7
alias python=python2.7

apt-get -y install python-setuptools python-dev build-essential 
apt-get -y install python-pip python-dev build-essential 
pip install pysam
apt-get -y install git-all

