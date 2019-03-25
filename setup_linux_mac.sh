#Get ViFi
git clone https://github.com/namphuon/ViFi.git
cd ViFi
VIFI_DIR=`pwd`

#Get data repos
wget https://raw.githubusercontent.com/circulosmeos/gdown.pl/master/gdown.pl
perl gdown.pl "https://drive.google.com/open?id=0ByYcg0axX7udUDRxcTdZZkg0X1k" data_repo.tar.gz
tar zxf data_repo.tar.gz
perl gdown.pl "https://drive.google.com/open?id=0Bzp6XgpBhhghSTNMd3RWS2VsVXM" data.zip
unzip data.zip

#Set up environmental variables
echo export VIFI_DIR=$VIFI_DIR >> ~/.bashrc
echo export AA_DATA_REPO=$PWD/data_repo >> ~/.bashrc
echo export REFERENCE_REPO=$PWD/data >> ~/.bashrc
source ~/.bashrc

#Set up reference for alignment
cat $AA_DATA_REPO//hg19/hg19full.fa $REFERENCE_REPO/hpv/hpv.unaligned.fas > $REFERENCE_REPO/hpv/hg19_hpv.fas

#Pull the Docker file
docker pull docker.io/namphuon/vifi 

docker run -v $REFERENCE_REPO/hpv/:/home/hpv/ docker.io/namphuon/vifi bwa index /home/hpv/hg19_hpv.fa

#Build reduced list of HMMs for testing
ls $VIFI_DIR/data/hpv/hmms/hmmbuild.[0-9].hmm > $VIFI_DIR/data/hpv/hmms/hmms.txt
source ~/.bashrc

#Run ViFi under docker mode on test dataset on reduced HMM list set
python $VIFI_DIR/scripts/run_vifi.py --cpus 2 --hmm_list $VIFI_DIR/data/hpv/hmms/hmms.txt -f $VIFI_DIR/test/data/test_R1.fq.gz -r $VIFI_DIR/test/data/test_R2.fq.gz -o $VIFI_DIR/tmp/docker/ --docker
