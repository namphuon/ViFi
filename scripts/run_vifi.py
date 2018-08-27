import sys, os, time
import argparse
import tempfile,shutil
import hg19util as hg19



def parse_args(reference_dir):
  parser = argparse.ArgumentParser()
  requiredNamed = parser.add_argument_group('Required named arguments')                                            
  requiredNamed.add_argument('-f', '--forward', type=str, action='store',
                      help='''forward FASTQ file''', required=True)
  requiredNamed.add_argument('-r', '--reverse', type=str, action='store', 
                      help='''reverse FASTQ file''', required=True)                                              
  outputOptions = parser.add_argument_group('Output options')                                              
  outputOptions.add_argument('-o', '--output_dir', default=".",
                      help='''output directory for all files created during run (default: current directory)''')
  outputOptions.add_argument('-p', '--prefix', default="output",
                      help='''prefix for output files''')                      
  runningOptions = parser.add_argument_group('Running options')                                              
  runningOptions.add_argument('-v', '--virus', default="hpv",
                      help='''VIRUS model to use.  (default: hpv)''')  
  runningOptions.add_argument('-c', '--cpus', default=1, type=int,
                      help='''number of cpus (default: 1)''')
  default_reference = "%s/<VIRUS>/hg19_<VIRUS>.fas" % (reference_dir)
  runningOptions.add_argument('--reference', default=default_reference, type=str,
                      help='''BWA reference (default: %s, VIRUS default is hpv)''' % default_reference)                                                                                   
  advOptions = parser.add_argument_group('Advanced options')                                              
  advOptions.add_argument('-C', '--chromosome_list',  dest='chromosome_list',  help='Optional file that contains list of chromosomes from reference organism.  All other sequences in the BAM file that are not in this list would be considered viral.  The file contains a single line with the space-delimited list of chromosomes belonging to the reference organism.  By default, all chromosomes starting with chr are considered human.',  metavar='FILE',  action='store',  type=str, default=None,)
  advOptions.add_argument('-b', '--bamfile', default=None,
                      help='''Use an existing NAME SORTED bamfile that might have been aligned with a different method.  If chromosome_list option not used, will assume reference genomes start with chr (default: None)''')
  advOptions.add_argument('-l', '--hmm_list', default=None,
                      help='''List of HMMs to include in search.  Useful if user wants to search against customized list of HMMs (default: None)''')                          
  advOptions.add_argument('-d', '--disable_hmms', default=False, type=bool,
                      help='''Disabled use of HMMs in search.  Useful if user only wants to use reference-based mapping to detect integrations.''')                                                
  options = parser.parse_args()
  return options

if __name__ == '__main__': 
  start_time = time.time()
  reference_dir = os.environ['REFERENCE_REPO']
  options = parse_args(reference_dir)  
  virus=options.virus
  if options.reference.find('<VIRUS>'):
     options.reference = options.reference.replace('<VIRUS>',virus)
  vifi_dir = os.environ['VIFI_DIR']
  if not os.path.exists(options.output_dir):
    os.makedirs(options.output_dir)  
  #Run BWA on input FASTQ files
  if options.bamfile is None:
    print "[Running BWA]: %f" % (time.time()-start_time)
    os.system("bwa mem -t %d -M %s %s %s | samtools view -bS - > %s/%s.bam" % (options.cpus, options.reference, options.forward, options.reverse, options.output_dir, options.prefix))
    print "[Finished BWA]: %f" % (time.time()-start_time)
    options.bamfile = "%s/%s.bam" % (options.output_dir, options.prefix)
    
  #Identify transitive reads
  os.chdir(options.output_dir)
  print "[Identifying chimeric reads]: %f" % (time.time()-start_time)  
  os.system("python %s/scripts/get_trans_new.py --unknown %s.unknown.bam --data %s --trans %s.trans.bam --viral %s.viral.bam %s" % (vifi_dir, options.prefix, options.bamfile, options.prefix, options.prefix, "" if options.chromosome_list is None else "--chrom %s" % options.chromosome_list))
  print "[Finished identifying chimeric reads]: %f" % (time.time()-start_time)  

  
  #Run HMMs, either use default HMMs found in directory or use list of HMMs given by user.  
  if options.disable_hmms == False:
    print "[Running HMMS]: %f" % (time.time()-start_time)
    if options.hmm_list is None:
      os.system("ls %s/%s/hmms/*.hmm > hmms.txt" % (reference_dir, virus))
      options.hmm_list = 'hmms.txt'
    if not os.path.exists('tmp/'):
      os.mkdir('tmp/')  
    os.system("python %s/scripts/run_hmms.py -t %d -b %s.unknown.bam -d tmp -H %s" % (vifi_dir, options.cpus, options.prefix, options.hmm_list))
    print "[Finished running HMMS]: %f" % (time.time()-start_time)
    
  #Cluster reads
  print "[Cluster and identify integration points]: %f" % (time.time()-start_time)
  os.system("python %s/scripts/merge_viral_reads.py --unknown %s.unknown.bam --trans %s.trans.bam --reduced tmp/temp/reduced.csv --map tmp/temp/unmapped.map --output %s.fixed.trans.bam" % (vifi_dir, options.prefix, options.prefix, options.prefix))
  os.system("samtools sort -m 2G -@ %d %s.fixed.trans.bam > %s.fixed.trans.cs.bam" % (options.cpus, options.prefix, options.prefix))
  os.system("samtools index %s.fixed.trans.cs.bam" % options.prefix)  
  os.system("samtools sort -m 2G -@ %d %s.viral.bam > %s.viral.cs.bam" % (options.cpus, options.prefix, options.prefix))  
  os.system("samtools index %s.viral.cs.bam" % options.prefix)
  os.system("python %s/scripts/cluster_trans_new.py --data %s.fixed.trans.cs.bam --output %s.clusters.txt %s" % (vifi_dir, options.prefix, options.prefix, "" if options.chromosome_list is None else "--chrom %s" % options.chromosome_list))  
  print "[Finished cluster and identify integration points]: %f" % (time.time()-start_time)
