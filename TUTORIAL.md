# ViFi Tutorial

We provide a brief tutorial on how to run ViFi and outline the different options below.  If there are further questions, please email Nam-phuong Nguyen (ndn006@eng.ucsd.edu) and we will update the tutorial with additional examples.  This tutorial assumes that Docker has been installed and ViFi has been downloaded.

[ViFi Structure](#ViFi Structure)
[Running ViFi to detect HPV integration from known strains](#Running-ViFi-to-detect-HPV-integration-from-known-strains)
[Running ViFi to detect HPV integration from novel strains](#Running-ViFi-to-detect-HPV-integration-from-novel-strains)
[Building reference package for different viral family](#Building-reference-package-for-different-viral-family)
[ViFi Output](#vifi-output)

### ViFi Structure

In this section, we outline the key components of ViFi.  First, ViFi requires a reference directory in order to function.  We point to the location of this directory using the $REFERENCE_REPO environmental variable.  Within this directory, the viral families of interest are separated out by directory.  For example, the default directory contains the following folders:

```
hpv
hbv
```

When one runs ViFi, the default virus is hpv, however, if the user specifies a different virus name (such as ebv), ViFi will assume that there is a folder with the same name within the $REFERENCE_REPO directory.  Within this folder, there several key files that ViFi expects to exist.  First, is a BWA index of the set of sequences from the human and viral references.  By default, the naming convention is hg19_<virus_name>.fas.  

### Running ViFi to detect HPV integration from novel strains

By default, ViFi will run using the hg19+HPV reference package for detecting integration using the HMMs.  Using the test data as an example, the user can then invoke:

```
python $VIFI_DIR/scripts/run_vifi.py -f $VIFI_DIR/test/test_R1.fq.gz -r $VIFI_DIR/test/test_R2.fq.gz --docker
```

The resulting integrations will be displayed in the folder in which this command was invoked in the file output.clusters.txt.  See [ViFi Output](#vifi-output) for more details on the output.  Reads that match to the viral HMMs will be labeled as viral_hmm within the output folder.

### Running ViFi to detect HPV integration from known strains

If the virus within the sample comes from a known strain (HPV16, for example), one can run ViFi without the HMMs.  To disable the use of HMMs within the search, the user can then invoke:

```
python $VIFI_DIR/scripts/run_vifi.py -f $VIFI_DIR/test/test_R1.fq.gz -r $VIFI_DIR/test/test_R2.fq.gz --docker --disable_hmms
```

The resulting integrations will be displayed in the folder in which this command was invoked in the file output.clusters.txt.  See [ViFi Output](#vifi-output) for more details on the output.


### ViFi Output

This will output the results of ViFi directly into the folder in which this command was invoked.  The output will be the following files (assuming that the prefix option is set to output):

output.bam                 output.fixed.trans.bam         output.trans.bam    output.viral.cs.bam
output.clusters.txt        output.fixed.trans.cs.bam      output.unknown.bam  output.viral.cs.bam.bai
output.clusters.txt.range  output.fixed.trans.cs.bam.bai  output.viral.bam

The output of ViFi is the list of read clusters discovered, and for each read cluster, the relaxed, stringent, and exact (if split reads are present) ranges are reported, aswell as the read names of the reads in the cluster.

The main output files of interest are
- output.clusters.txt
- output.clusters.txt.range

output.clusters.txt is a tab delimited file that reports the human integration range, the number of reads supporting the integration, and the number of reads mapped to the forward/reverse strand of the human region, as well as the number of viral reads mapping to the virus sequence.  It also includes the names of each discordant read supporting the integration.

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

output.clusters.txt.range is a much more condensed summary of the results, showing just the integration range on the
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

Finally, ViFi outputs several working files that can be deleted after a run.  These are:
1. hmms.txt - The list of HMM files used during the run
2. output.bam - The aligned (name-sorted order) BAM file containing the input reads
3. output.unknown.bam - A BAM file containing all paired reads in which one or both paired end reads that did not align to any known reference.  ViFi will then search these reads against the HMMs to identify any viral reads.
4. output.viral.bam - A BAM file containing all paired reads that only aligned to viral references
5. output.viral.cs.bam - A coordinate sorted BAM file containing all paired reads that only aligned to viral references
6. output.trans.bam - A BAM file containing all paired reads in which one read aligned to the human and the other aligned to the viral reference.
7. output.fixed.trans.bam - A BAM file created by merging 6. and any human/viral paired end reads discovered by running the viral HMMs on 3.
8. output.fixed.trans.cs.bam - A coordinate sorted BAM file of 7.


