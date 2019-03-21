import sys, os, time
import argparse
import tempfile,shutil
import hg19util as hg19



def parse_args(reference_dir):
  parser = argparse.ArgumentParser()
  requiredNamed = parser.add_argument_group('Required named arguments')                                              
  requiredNamed.add_argument('-f', '--forward', type=str, action='store',
                      help='''forward FASTQ file''')
  requiredNamed.add_argument('-r', '--reverse', type=str, action='store', 
                      help='''reverse FASTQ file''')     
  requiredNamed.add_argument('-b', '--bamfile', default=None,
                      help='''Use an existing NAME SORTED bamfile that might have been aligned with a different method.  If chromosome_list option not used, will assume reference genomes start with chr (default: None)''')
                                                               
  outputOptions = parser.add_argument_group('Output options')                                              
  outputOptions.add_argument('-o', '--output_dir', default=".",
                      help='''output directory for all files created during run (default: current directory)''')
  outputOptions.add_argument('-p', '--prefix', default="output",
                      help='''prefix for output files''')                      
  runningOptions = parser.add_argument_group('Running options')                                              
  runningOptions.add_argument('-v', '--virus', default="hpv",
                      help='''VIRUS model to use.  (default: hpv)''')                        
  runningOptions = parser.add_argument_group('Running options')                                              
  runningOptions.add_argument('--docker', action='store_true',
                      help='''Run via docker''')
  runningOptions.add_argument('-c', '--cpus', default=1, type=int,
                      help='''number of cpus (default: 1)''')
  default_reference = "%s/<VIRUS>/hg19_<VIRUS>.fas" % (reference_dir)
  runningOptions.add_argument('--reference', default=default_reference, type=str,
                      help='''BWA reference (default: %s, VIRUS default is hpv)''' % default_reference)                                                                                   
  advOptions = parser.add_argument_group('Advanced options')                                              
  advOptions.add_argument('-C', '--chromosome_list',  dest='chromosome_list',  help='Optional file that contains list of chromosomes from reference organism.  All other sequences in the BAM file that are not in this list would be considered viral.  The file contains a single line with the space-delimited list of chromosomes belonging to the reference organism.  By default, the 23 full-length chromosomes from hg19 are considered human.',  metavar='FILE',  action='store',  type=str, default=None,)
  advOptions.add_argument('-l', '--hmm_list', default=None,
                      help='''List of HMMs to include in search, one HMM filename per line.  If running under Docker mode, all the hmms must be within the same folder as the hmm_list file.  Useful if user wants to search against customized list of HMMs (default: None)''')                          
  advOptions.add_argument('-d', '--disable_hmms', default=False, type=bool,
                      help='''Disabled use of HMMs in search.  Useful if user only wants to use reference-based mapping to detect integrations.''')                                         
  options = parser.parse_args()
  
  #Build input_string for docker run as well as vifi_args, as well as verify all options
  input_string = ""
  vifi_string = ""
  
  #This block verifies the existence of the input data, as well as making sure that either
  #the BAM file or the FASTQ files are given as input
  if options.bamfile is not None and (options.forward is not None or options.reverse is not None) :     
    parser.error('Either supply a BAM file (-b) or both fastq files (-r/-f)')
  elif options.bamfile is None and (options.forward is None or options.reverse is None):
    parser.error('Must supply a BAM file (-b) or both fastq files (-r/-f)')
  elif options.bamfile is None and (options.forward is not None and options.reverse is not None):
    if not os.path.exists(options.forward) or not os.path.exists(options.reverse):
      parser.error('Unable to find %s or %s' % (options.forward, options.reverse))
    else:
      forward_dir = os.path.dirname(os.path.realpath(options.forward))
      reverse_dir = os.path.dirname(os.path.realpath(options.reverse))
      input_string+= "-v %s:/home/fastq1/ -v %s:/home/fastq2/ -e READ1=%s -e READ2=%s " % (forward_dir, reverse_dir, os.path.basename(options.forward), os.path.basename(options.reverse))    
      vifi_string+= "-f /home/fastq1/%s -r /home/fastq2/%s " % (os.path.basename(options.forward), os.path.basename(options.reverse))
  elif options.bamfile is not None and (options.forward is None or options.reverse is None):
    if not os.path.exists(options.bamfile):
      parser.error('Unable to find %s' % (options.bamfile))
    else:
      bam_dir = os.path.dirname(os.path.realpath(options.bamfile))
      input_string+= "-v %s:%s-e BAMFILE=%s " % (bam_dir, '/home/bam/', os.path.basename(options.bamfile))                
      vifi_string+= "-b /home/bam/%s " % (os.path.basename(options.bamfile))
  if options.chromosome_list is not None:
    if not os.path.exists(options.chromosome_list):
      parser.error('Unable to find %s' % (options.chromosome_list))
    else:
      chromosome_dir = os.path.dirname(os.path.realpath(options.chromosome_list))
      chromosome = os.path.basename(options.chromosome_list) 
      input_string+= "-v %s:%s -e CHROMOSOME=%s " % (chromosome_dir, '/home/chromosomes/', os.path.basename(options.chromosome))
      vifi_string+= "-C /home/chromosomes/%s " % (os.path.basename(options.chromosome))                
  if options.hmm_list is not None:
    if not os.path.exists(options.hmm_list):
      parser.error('Unable to find %s' % (options.hmm_list))
    else:
      hmm_list_dir = os.path.dirname(os.path.realpath(options.hmm_list))
      hmm_list = os.path.basename(options.hmm_list)
      if options.docker == True:     
        temp_list = create_new_hmm_list(hmm_list_dir, hmm_list)
        input_string+= "-v %s:%s -e HMM_LIST=%s " % (hmm_list_dir, '/home/hmm_list/', os.path.basename(temp_list.name))                
        vifi_string+= "-l /home/hmm_list/%s " % (os.path.basename(temp_list.name))      
        options.temp_list = temp_list.name          
  if options.reference.find('<VIRUS>') == -1:
    if not os.path.exists(options.default_reference):
      parser.error('Unable to find %s' (options.default_reference))
    else:
      #Fix default reference
      vifi_string+= "-reference %s " % (options.default_reference)                   
  input_string+= "-v %s:/home/repo/data " % os.environ['REFERENCE_REPO']
  input_string+= "-v %s:/home/output/  -e REFERENCE_REPO=/home/repo/data "  % (os.path.realpath(options.output_dir))
  input_string+= "-v %s:/home/data_repo/ " % os.environ['AA_DATA_REPO']
  vifi_string+= "-c %d " % (options.cpus)               
  vifi_string+= "-o /home/output/ -p %s " % (options.prefix)               
  options.cmd_string = 'docker run %s docker.io/namphuon/vifi python "scripts/run_vifi.py" %s' % (input_string, vifi_string)  
  #"docker run -e CPUS=$CPUS -v $REFERENCE_REPO:/home/repo/data -v $INPUT_DIR:/home/fastq/ -e READ1=$READ1 -e READ2=$READ2 -v $AA_DATA_REPO:/home/data_repo/ -v $OUTPUT_DIR:/home/output/ vifi:latest python "scripts/run_vifi.py %s" % (input_string  
  return options

def create_new_hmm_list(hmm_list_dir, hmm_list):
  input = open("%s/%s" % (hmm_list_dir, hmm_list), 'r')
  output = open("%s/%s" % (hmm_list_dir, next(tempfile._get_candidate_names())), 'w')
  for line in input:
    line = line.strip()
    if os.path.exists(line):
      name = os.path.basename(line)
      output.write("%s/%s\n" % ("/home/hmm_list/", name))
  output.close()    
  os.system('chmod +0777 %s' % output.name)
  return output

def run_docker(options):  
  os.system(options.cmd_string)
  os.system('rm %s' % options.temp_list)
  
if __name__ == '__main__': 
  start_time = time.time()
  reference_dir = os.environ['REFERENCE_REPO']
  options = parse_args(reference_dir)
  if options.docker == True:
    run_docker(options)
    exit(0)
  virus=options.virus
  if options.reference.find('<VIRUS>'):
     options.reference = options.reference.replace('<VIRUS>',virus)
  vifi_dir = os.environ['VIFI_DIR']
  if not os.path.exists(options.output_dir):
    os.makedirs(options.output_dir)
  options.output_dir = os.path.abspath(options.output_dir)
  #Run BWA on input FASTQ files
  if options.bamfile is None:
    print "[Running BWA]: %f" % (time.time()-start_time)
    os.system("bwa mem -t %d -M %s %s %s | samtools view -bS - > %s/%s.bam" % (options.cpus, options.reference, options.forward, options.reverse, options.output_dir, options.prefix))
    print "[Finished BWA]: %f" % (time.time()-start_time)
    options.bamfile = "%s/%s.bam" % (options.output_dir, options.prefix)
  if options.chromosome_list is not None:
    options.chromosome_list = os.path.abspath(options.chromosome_list)
  if options.hmm_list is not None:
    options.hmm_list = os.path.abspath(options.hmm_list)
  
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
