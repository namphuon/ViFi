#!/usr/bin/python
# -*- coding: utf-8 -*-

import pysam
import argparse
from time import clock
from collections import Counter

from sets import Set

def is_float(value):
  try: 
    float(value)
    return True
  except ValueError:
    return False
    
def read_scores_file(hmm_file):
  input = open(hmm_file, 'r')  
  scores = dict()
  for line in input:
    results = line.strip().split(',')
    scores[results[0]] = [float(x) if is_float(x) else x for x in results]
  input.close()
  return scores

def read_map(map_file, scores, delimiter = '\t', direction = 'reverse'):
  input = open(map_file, 'r')  
  map = dict()
  for line in input:
    results = line.strip().split(delimiter)    
    if results[1] in scores:      
      results[0] = results[0].replace('/1','').replace('/2','').replace('@','')
      if results[0] in map and map[results[0]][2] < scores[results[1]][2]:
        map[results[0]] = scores[results[1]]
      elif results[0] not in map:
        map[results[0]] = scores[results[1]]
  input.close()
  return map


parser = \
  argparse.ArgumentParser(description='Merges viral reads identified via HMMs'
              )

parser.add_argument(
  '--trans',
  dest='transName',
  help='Output BAM file of trans reads between hg19 and viral',
  metavar='FILE',
  action='store',
  type=str,
  nargs=1,
  )
  

parser.add_argument(
  '--unknown',
  dest='unknownName',
  help='BAM file of reads with only one read mapped to between hg19 and other to unknown',
  metavar='FILE',
  action='store',
  type=str,
  nargs=1,
  )

parser.add_argument(
  '--threshold',
  dest='threshold',
  help='E-value threshold to accept as virus read',
  metavar='float',
  action='store',
  type=float,
  default=0.01,
  )

  
parser.add_argument(
  '--reduced',
  dest='reducedName',
  help='CSV file containing HMM scores',
  metavar='FILE',
  action='store',
  type=str,
  nargs=1,
  )

parser.add_argument(
  '--map',
  dest='mapName',
  help='Map file containing mapping of scores and sequences',
  metavar='FILE',
  action='store',
  type=str,
  nargs=1,
  )
  
parser.add_argument(
  '--output',
  dest='outputName',
  help='Output BAM file of updated trans reads',
  metavar='FILE',
  action='store',
  type=str,
  nargs=1,
  )
args = parser.parse_args()

transFile = pysam.Samfile(args.transName[0], 'rb')
unknownFile = pysam.Samfile(args.unknownName[0], 'rb')

references = transFile.header
references['SQ'].append({'LN': 5000, 'SN': 'viral_hmm'})

outputFile = pysam.Samfile(args.outputName[0], 'wb', header=references)
scores = read_scores_file(args.reducedName[0])
mapping = read_map(args.mapName[0], scores)  
outputFile.header['SQ'] = references
for read in unknownFile.fetch(until_eof=True):
  if read.qname in mapping and mapping[read.qname][1] < args.threshold:
    if read.is_unmapped:      
      read.tid = len(references['SQ'])-1
      read.is_unmapped = False
      read.pos = 1      
      read.cigartuples = [(0,read.qlen)]
#       matched = int(abs(mapping[read.qname][6]-mapping[read.qname][5]))+1
#       if matched == read.qlen:
#         read.cigartuples = [(0,matched)]
#       else:
#         read.cigartuples = [(0,matched),(4,read.qlen-matched)]
    else:
      read.mate_is_unmapped = False
      read.mpos = 1      
      read.mrnm = len(references['SQ'])-1
    outputFile.write(read)

for read in transFile.fetch(until_eof=True):
  outputFile.write(read)

outputFile.close()

