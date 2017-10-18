import sys, os, time
import argparse
import tempfile,shutil
import hg19util as hg19



def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('-v', '--virus', default="hpv",
                      help='''virus model (default: hpv)''')
  parser.add_argument('-f', '--forward', default="",
                      help='''forward FASTQ file''')
  parser.add_argument('-r', '--reverse', default="",
                      help='''reverse FASTQ file''')                      
  parser.add_argument('-o', '--output_dir', default=".",
                      help='''output directory for all files created during run (default: current directory)''')
  parser.add_argument('-p', '--prefix', default="output",
                      help='''prefix for output files''')                      
  parser.add_argument('-c', '--cpus', default=1,
                      help='''number of cpus (default: 1)''')                                                              
  options = parser.parse_args()
  return options
    
    
if __name__ == '__main__': 
  start_time = time.time()
  options = parse_args()
  reference_dir = os.environ['REFERENCE_REPO']
  vifi_dir = os.environ['VIFI_DIR']
  if !os.path.exists(options.output_dir):
    os.makedirs(options.output_dir)
  
  #Run BWA on input FASTQ files
  print "[Running BWA]: %f" % time.time()-start_time
  os.system("bwa mem -t %d -M %s/hg19_%s.fas %s %s | samtools view -bS - > %s/hg19_%s.bam" % (options.cpus, options.output_dir, options.virus, options.forward, options.reverse, options.output_dir, options.virus))
  print "[Finished BWA]: %f" % time.time()-start_time
  
  
  #Identify transitive reads
  os.chdir(options.output_dir)
  print "[Identifying chimeric reads]: %f" % time.time()-start_time  
  os.system("python %s/scripts/get_trans_new.py --unknown hg19_%s.unknown.bam --data hg19_%s.bam --trans hg19_%s.trans.bam --viral hg19_%s.viral.bam" % (vifi_dir, virus, virus, virus, virus))
  print "[Finished identifying chimeric reads]: %f" % time.time()-start_time  

  #Run HMMs
  print "[Running HMMS]: %f" % time.time()-start_time
  os.system("ls %s/%s/hmms/*.hmm > hmms.txt" % (reference_dir, virus)) 
  os.system("perl -I %s/lib %s/scripts/identify_viruses.pl -b hg19_%s.unknown.bam -d -H hmms.txt" % (vifi_dir, virus))
  print "[Finished running HMMS]: %f" % time.time()-start_time
  
  #Cluster reads
  print "[Cluster and identify integration points]: %f" % time.time()-start_time
  os.system("python %s/scripts/merge_viral_reads.py --unknown hg19_%s.unknown.bam --trans hg19_%s.trans.bam --reduced test/temp/reduced.csv --map test/temp/unmapped.map --output hg19_%s.fixed.trans.bam" % (vifi_dir, virus, virus, virus))
  os.system("samtools sort -m 2G -@ %d hg19_%s.fixed.trans.bam > hg19_%s.fixed.trans.cs.bam" % (options.cpus, virus, virus))
  os.system("samtools sort -m 2G -@ %d hg19_%s.viral.bam > hg19_%s.viral.cs.bam" % (options.cpus, virus, virus))
  os.system("samtools index hg19_%s.fixed.trans.cs.bam" % virus)
  os.system("samtools index hg19_%s.viral.cs.bam" % virus)
  os.system("python %s/scripts/cluster_trans_new.py --data hg19_%s.fixed.trans.cs.bam --output %s.clusters.txt" % (vifi_dir, virus, prefix))  
  print "[Finished cluster and identify integration points]: %f" % time.time()-start_time