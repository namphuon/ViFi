import sys, os, time,pysam
import argparse
import tempfile,shutil

def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('-t', '--threads', type=int, action='store',
                      help='''Threads to use''', required=False, default = 1)
  parser.add_argument('-b', '--bamfile', type=str, action='store', 
                      help='''Bamfile to analyze''', required=True)                                              
  parser.add_argument('-d', '--directory', default=".",
                      help='''output directory for all files created during run (default: current directory)''')
  parser.add_argument('-H', '--hmm_list', default=None, type=str, required=True,
                      help='''File containing HMMs to use''')                      
  options = parser.parse_args()
  return options

def read_hmm_file(hmm_list):
  hmms = {}
  input = open(hmm_list, "r");
  idx = 0;
  for line in input:
    #Todo, instead of index, note which viral family the hmm came from
    hmms[idx] = line.strip()
    idx+=1
  input.close()
  return hmms

def run_pipeline(options):
  hmms = read_hmm_file(options.hmm_list)
  if not os.path.exists('%s/logs' % options.directory):
    os.mkdir('%s/logs' % options.directory)
  if not os.path.exists('%s/temp' % options.directory):
    os.mkdir('%s/temp' % options.directory)
  if not os.path.exists('%s/temp/unmapped.fas'):
    prepare_unmapped_sequences(options)
  #Now search HMMs against reads
  print "Running HMMs"
  start_time = time.time() 
  for i in hmms.keys():
    print "\tRunning HMM %s" % hmms[i]
    os.system('nhmmer -o %s/temp/hmmsearch.%s --noali --cpu %d %s %s/temp/unmapped.fas' % (options.directory, str(i), options.threads, hmms[i], options.directory))
  end_time = time.time() - start_time;
  print "Finished running against HMMs: %fs" % end_time    
  print "Processing results\n";
  scores = {}
  for i in hmms.keys():
    scores = read_nhmmer_result("%s/temp/hmmsearch.%s" % (options.directory, i), scores)
  output = open('%s/temp/reduced.csv' % options.directory, 'w')
  for score in scores.keys():
    output.write("%s\n" % (','.join([str(s) for s in scores[score]])))
  output.close()

def read_nhmmer_result(file, scores):
  input = open(file, 'r')
  start_line = 'Query:'
  start = False
  for line in input:
    if start == False and line.find(start_line) == 0:
      start = True
      foo = input.next()
      foo = input.next()
      foo = input.next()      
    if start == True and line.find('read_') != -1:
      res = line.split()
      if scores.setdefault(res[3],[res[3], res[0], float(res[1])])[2] < float(res[1]):
        scores[res[3]] = [res[3], res[0], float(res[1])]
    elif line.strip() == "" and start == True:
      break
  return scores

def prepare_unmapped_sequences(options):
  start_time = time.time() 
  counter = 0;
  bam = pysam.Samfile(options.bamfile, 'rb')
  fas = open('%s/temp/unmapped.fas' % options.directory, 'wb')
  map = open('%s/temp/unmapped.map' % options.directory, 'wb')
  
  for read in bam:
    if read.is_unmapped and not read.mate_is_unmapped and bam.references[read.rnext].find('chr') == 0:
      fas.write('>read_%d\n%s\n' % (counter, read.seq))
      map.write('%s\tread_%d\n' % (read.qname, counter))
      counter+=1
  fas.close()
  map.close()
  bam.close()
  end_time = time.time() - start_time;
  print "Prepared sequences for searching against HMMs: %fs" % end_time

if __name__ == '__main__': 
  start_time = time.time()
  reference_dir = os.environ['REFERENCE_REPO']
  options = parse_args()  
#   class Object(object):
#       pass
#   options = Object()  
#   options.bamfile='output.unknown.bam'
#   options.directory='tmp'
#   options.hmm_list='hmms.txt'
#   options.threads=1
  run_pipeline(options)