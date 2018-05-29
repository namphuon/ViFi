# ViFi

ViFi is a tool for detecting viral integration and fusion mRNA sequences from Next Generation Sequencing data.  Unlike standard approaches that use reference-based read mapping for identification of viral reads, ViFi uses both reference-based read mapping and a phylogenetic-based approach to identify viral reads.  ViFi also incorporates mappability scores of the reads to filter out false positive integration detection.  The end result is a tool that can accurately and precisely detect integrated viruses, even if the viruses are highly mutated or novel strains.

ViFi is currently in alpha testing, is is constantly undergoing revisions.  High on the priority list is an easier installation process, as well as improve user interface.  Please report any problems/bugs to Nam Nguyen (ndn006@eng.ucsd.edu) so that ViFi can be improved and problems can be quickly corrected.  


## Installation:
We provide instructions for installing ViFi on Linux below.  

1. ViFi download (if you have not already cloned this source code):
```
git clone https://github.com/namphuon/ViFi.git
```
2. Install Dependencies:
   1. Python 2.7
   ```
   sudo dnf install python2
   ```
   2. Pysam verion 0.9.0 or higher (https://github.com/pysam-developers/pysam):
   ```
   sudo pip install pysam
   ```
   3. Samtools 1.3.1 or higher (www.htslib.org/)
   ```
   sudo apt-get install samtools
   ```
   4. BWA 0.7.15 or higher (bio-bwa.sourceforge.net/)
   ```
   sudo apt-get install bwa
   ```
   5. Install HMMER v3.1b2 and have it on the path (http://hmmer.org/)
   ```
   sudo apt-get install hmmer
   ```
3. Set the ViFi directory and include the python source to your Python path
```
echo export VIFI_DIR=/path/to/ViFi >> ~/.bashrc
echo export PYTHONPATH=/path/to/ViFi:$PYTHONPATH >> ~/.bashrc
```
4. Download the data repositories:
While we include some annotations, we are unable to host some large files in the git repository.  These may be downloaded from https://drive.google.com/open?id=0ByYcg0axX7udUDRxcTdZZkg0X1k. Thanks to Peter Ulz for noticing incorrect link earlier.
```
tar zxf data_repo.tar.gz
echo export AA_DATA_REPO=$PWD/data_repo >> ~/.bashrc
source ~/.bashrc
```
5. Download the HMM models:
We have pre-build HMM models for HPV and HBV.  They can be downloaded from https://drive.google.com/open?id=0Bzp6XgpBhhghSTNMd3RWS2VsVXM.
```
unzip data.zip
echo export REFERENCE_REPO=$PWD/data/ >> ~/.bashrc
```
6.  Build a BWA index on the reference sequences from human+viral sequences:
We show an example of building an index of human+viral sequences using Hg19 and **HPV** and **HBV** below.  However
any reference organism+viral family could be used.
```
cat $AA_DATA_REPO//hg19/hg19full.fa $REFERENCE_REPO/hpv/hpv.unaligned.fas > $REFERENCE_REPO/hpv/hg19_hpv.fas
bwa index $REFERENCE_REPO/hpv/hg19_hpv.fas

cat $AA_DATA_REPO//hg19/hg19full.fa $REFERENCE_REPO/hbv/hbv.unaligned.fas > $REFERENCE_REPO/hbv/hg19_hbv.fas
bwa index $REFERENCE_REPO/hbv/hg19_hbv.fas
```
## Running ViFi 

We show the most basic example of running ViFi below.  This version assumes that the user has
followed all the previous steps.  More advanced options, such as using a customized reference organism/viral
family is provided in the [Advanced Notes][#advanced_notes] section.  
```
python run_vifi.py -f <input_R1.fq.gz> -r <input_R2.fq.gz> -o <output_dir>
```

Note that this version defaults to searching for **HPV**.  To search for HBV, run the following command.
```
python run_vifi.py -f <input_R1.fq.gz> -r <input_R2.fq.gz> -o <output_dir> -v hbv
```

## ViFi Output

The output of ViFi is the list of read clusters discovered, and for each read cluster, the relaxed, stringent, and exact (if split reads are present) ranges are reported, aswell as the read names of the reads in the cluster.

The main output files of interest are
- \<prefix\>.clusters.txt
- \<prefix\>.clusters.txt.range

\<prefix\>.clusters.txt is a tab delimited file that reports the human integration range, the number of reads supporting the integration, and the number of reads mapped to the forward/reverse strand of the human region, as well as the number of viral reads mapping to the virus sequence.  It also includes the names of each discordant read supporting the integration.

Below is the sample:
```
#chr    minpos  maxpos  #reads  #forward        #reverse
##================================================================
chr19   36212224        36212932        7       4       3
##ERR093797.9977893     chr19   36212224        True    False
##ERR093797.7073606     chr19   36212403        True    True
```

The first line is the header information.  Afterward, each integration cluster is separated by a line
containing **=**.  The first line of an integration cluster describes the following:

1. Reference chromosome (chr19)
2. Minimum reference position of all mapped reads belonging to that cluster (36212224)
3. Maximum reference positions of all mapped reads belonging to that cluster (36212932)
4. Number of read pairs belonging to this cluster (7)
5. Number of reads mapped to the forward reference strand (4)
6. Number of reads mapped to the forward reference strand (3)

After this line, each read pair that mapped to this cluster is displayed.  The information is
1. Read name (ERR093797.9977893)
2. Reference chromosome (chr19)
3. Starting read map location (36212224)
4. Read is on the reverse strand (True)
5. Read is read1 (False)

\<prefix\>.clusters.txt.range is a much more condensed summary of the results, showing just the integration range on the
human reference (based upon discordant reads) and attempts to identify the exact integration point if split reads are available.

Below is a sample:
```
Chr,Min,Max,Split1,Split2
chr19,36212564,36212564,-1,-1
```

The first line is header information.  Afterward, each line is information about the cluster.  For example,

1. Reference chromosome (chr19)
2. Minimum reference position of all mapped reads belonging to that cluster (36212224)
3. Maximum reference positions of all mapped reads belonging to that cluster (36212932)
4. If split read exists, minimum split read mapped range, -1 if no split read exists (-1)
5. If split read exists, maximum split read mapped range, -1 if no split read exists (-1)

## Dockerized ViFi

We have also created a dockerized version of ViFi to enable easier time running.  The docker version of ViFi can be obtained
by installing Docker (https://www.docker.com/), and running the following command:

docker pull namphuon/vifi

To run the dockerized version of ViFi, first create the data repositories as above, including setting the environmental variables. 
Next, run the following script in the ViFi scripts directory:

`docker_vifi.sh <INPUT_DIR> <READ1>  <READ2>  <OUTPUT>  <CPUS>` 

where <INPUT_DIR> is the directory containing the <READ1> and <READ2> files, and <CPUS> is the number of
CPUs to use.  Note that the full path must
be given for the input and output directory, and the $AA_DATA_REPO and $REFERENCE_REPO variables must be
set in order for the script to find the necessary files.  

Example:

If /home/input/ contains read1.fastq.gz and read2.fastq.gz, then

sh docker_vifi.sh /home/input read1.fastq.gz read2.fastq.gz /home/output/ 2

## References
1. Nguyen ND, Deshpande V, Luebeck J, Mischel PS, Bafna V (2018) ViFi: accurate detection of viral integration and mRNA fusion reveals indiscriminate and unregulated transcription in proximal genomic regions in cervical cancer. Nucleic Acids Res (April):1–17.

# [Advanced Notes][#advanced_notes]

## Building evolutionary models

ViFi can be run with and without evolutionary models (i.e., the HMMs).  The HMMs

## Building Alignment and Tree on viral family of interest

ViFi can build HMMs from any viral family if there is an existing FASTA alignment and NEWICK tree on the viral
sequences.  Note that the sequences should be phylogenetically related to each other (i.e., do not mix HPV and HBV
sequences).  Any standard alignment method and tree reconstruction method can be used.  In our paper, we used [PASTA](https://github.com/smirarab/pasta) to construct our alignment and tree and provide the steps in doing this below.
Instructions on installing and running PASTA can be found [here](https://github.com/smirarab/pasta).

## Building HMMs

We created script to allow easy creation of the HMMs used within ViFi for a viral family of interesting.  To

Requires:
## 1) Python 2.7
## 2) Dendropy verion 4.0.0 or higher (https://github.com/jeetsukumaran/DendroPy):
sudo pip install dendropy

## Using customized reference 

If you want to use a customized reference or a reference for a different organism, you can inform ViFi 
of the reference sequences by supplying a chromosome file to ViFi using the **--chromosome_list**.  The file 
format is a single line that has the sequence names delimited by spaces.  For example:

```
mouse_chr1 mouse_chr2
```

would inform ViFi that any other sequences found in the BAM file that does not match mouse_chr1 and mouse_chr2 are
considered viral sequences.