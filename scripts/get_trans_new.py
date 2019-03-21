#!/usr/bin/python
# -*- coding: utf-8 -*-

import pysam
import argparse
from time import clock
from collections import Counter

from sets import Set

parser = \
  argparse.ArgumentParser(description='Extract reads with one end mapping to hg19 reference and other end mapping to viral genome'
              )

parser.add_argument(
  '--unknown',
  dest='unknownName',
  help='Namesorted unknown BAM file containing reads that were partially unalignable',
  metavar='FILE',
  action='store',
  type=str,
  nargs=1,
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
  '--data',
  dest='dataName',
  help='Input BAM file',
  metavar='FILE',
  action='store',
  type=str,
  nargs=1,
  )
  
parser.add_argument(
  '--viral',
  dest='viralName',
  help='Output BAM file of viral reads',
  metavar='FILE',
  action='store',
  type=str,
  nargs=1,
  )

parser.add_argument(
  '--misc',
  dest='miscName',
  help='Output BAM file of misc. double unmapped reads',
  metavar='FILE',
  action='store',
  default = None,
  type=str,
  nargs=1,
  )
    
parser.add_argument(
  '--output',
  dest='outputName',
  help='Output BAM file of augmented trans reads between reference organism and viral sequences',
  metavar='FILE',
  action='store',
  type=str,
  nargs=1,
  )

parser.add_argument(
  '--chrom_list',
  dest='chrom_list',
  help='Optional file that contains list of chromosomes from reference organism.  All other sequences in the BAM file that are not in this list would be considered viral.  The file contains a single line with the space-delimited list of chromosomes belonging to the reference organism.  By default, all chromosomes starting with chr are considered human.',
  metavar='FILE',
  action='store',
  type=str,
  default = None,
  nargs=1,
  )
  
opts = parser.parse_args()

bamFile = pysam.Samfile(opts.dataName[0], 'rb')
transFile = pysam.Samfile(opts.transName[0], 'wb', template=bamFile)
viralFile = pysam.Samfile(opts.viralName[0], 'wb', template=bamFile)
unknownFile = pysam.Samfile(opts.unknownName[0], 'wb', template=bamFile)
chrom_list = {}
foo = [chrom_list.setdefault('human' if ref.find('chr') == 0 else 'viral',Set()).add(ref) for ref in bamFile.references]
if opts.chrom_list is not None:
  chrom_list = {}
  input = open(opts.chrom_list[0], 'r')
  [chrom_list.setdefault('human',Set()).add(l) for l in input.next().strip().split()]
  #[chrom_list.setdefault('viral',Set()).add(l) for l in input.next().strip().split(' ')]
miscFile = None
if (opts.miscName is not None):
  miscFile = pysam.Samfile(opts.miscName[0], 'wb', template=bamFile)  

qname = ''
q1ref = Set([])
q2ref = Set([])
q1aligns = []
q2aligns = []

#hg19refs = Set(map(lambda x: 'chr' + str(x), range(1, 23) + ['X', 'Y',
#         'M']))

transReads = 0
totalReads = 0
viralReads = 0
unknownReads = 0
oldalign = ''
totalaligns = 0

for a in bamFile:
  totalaligns += 1
  if a.qname != qname and qname != '':
    if len([b for b in q1aligns if not b.is_secondary and not b.is_supplementary]) != 1 \
      or len([b for b in q2aligns if not b.is_secondary and not b.is_supplementary]) != 1:
      print 'No primary alignment', a.qname
      print [(b, str(b)) for b in q1aligns if not b.is_secondary]
      print [(b, str(b)) for b in q2aligns if not b.is_secondary]
      exit()
    trans = False
    viral = False
    unknown = False
    totalReads += 1
    if totalReads % 100000 == 0:
      print clock(), totalReads, 'reads done: #(Trans reads) =', \
        transReads, viralReads, a.qname, qname
    len_q1ref = len(chrom_list['human'].intersection(q1ref))
    len_q2ref = len(chrom_list['human'].intersection(q2ref))
    if len_q1ref > 0 \
      and len_q2ref == 0 and len(q2ref) \
      > 0:
      trans = True
    if len_q2ref > 0 \
      and len_q1ref == 0 and len(q1ref) \
      > 0:
      trans = True
    if len_q2ref == 0 \
      and len_q1ref == 0 and (len(q1ref)
        > 0 or len(q2ref) > 0):
      viral = True
    if trans == True:
      transReads += 1
      if len([b for b in q1aligns if bamFile.getrname(b.tid)
           in chrom_list['human']]) <= 3 or len([b for b in q2aligns
          if bamFile.getrname(b.tid) in chrom_list['human']]) <= 3:
        for b in q1aligns + q2aligns:
          transFile.write(b)
    if viral == True:
      viralReads += 1
      for b in q1aligns + q2aligns:
        viralFile.write(b)
    if viral == False and trans == False:
      #At least one read maps to human and at least one read is unmapped
      if (len([read for read in q1aligns if not read.is_unmapped and bamFile.references[read.tid] in chrom_list['human']])+len([read for read in q2aligns if not read.is_unmapped and bamFile.references[read.tid] in chrom_list['human']])) > 0 and (len([read for read in q2aligns if read.is_unmapped]) + len([read for read in q1aligns if read.is_unmapped])) > 0:
        unknownReads += 1
        for b in q1aligns + q2aligns:
          seq = Counter(q1aligns[0].seq + q2aligns[0].seq)
          if 'N' in seq and seq['N'] >= 5:
            continue
          unknownFile.write(b)
      if (miscFile is not None and len([read for read in q2aligns if read.is_unmapped]) + len([read for read in q1aligns if read.is_unmapped])) >= 2:
        for b in q1aligns + q2aligns:
          seq = Counter(q1aligns[0].seq + q2aligns[0].seq)
          if 'N' in seq and seq['N'] >= 5:
            continue
          miscFile.write(b)
      

  if a.qname != qname:
    # print qname, hg19refs, q1ref, q2ref, hg19refs.intersection(q1ref), hg19refs.intersection(q2ref)
    qname = a.qname
    q1aligns = []
    q2aligns = []
    q1ref = Set([])
    q2ref = Set([])
  if a.is_read1:
    q1aligns.append(a)
    if not a.is_unmapped and not a.mate_is_unmapped and a.tid != -1:
      q1ref.add(bamFile.getrname(a.tid))
  else:
    q2aligns.append(a)
    if not a.is_unmapped and not a.mate_is_unmapped and a.tid != -1:
      q2ref.add(bamFile.getrname(a.tid))
  oldalign = a

trans = False
viral = False
totalReads += 1

len_q1ref = len(chrom_list['human'].intersection(q1ref))
len_q2ref = len(chrom_list['human'].intersection(q2ref))
if len_q1ref > 0 \
  and len_q2ref == 0 and len(q2ref) \
  > 0:
  trans = True
if len_q2ref > 0 \
  and len_q1ref == 0 and len(q1ref) \
  > 0:
  trans = True
if len_q2ref == 0 \
  and len_q1ref == 0 and (len(q1ref)
    > 0 or len(q2ref) > 0):
  viral = True
if trans == True:
  transReads += 1
  if len([b for b in q1aligns if bamFile.getrname(b.tid)
       in chrom_list['human']]) <= 3 or len([b for b in q2aligns
      if bamFile.getrname(b.tid) in chrom_list['human']]) <= 3:
    for b in q1aligns + q2aligns:
      transFile.write(b)
if viral == True:
  viralReads += 1
  for b in q1aligns + q2aligns:
    viralFile.write(b)
if viral == False and trans == False:
  #At least one read maps to human and at least one read is unmapped
  if (len([read for read in q1aligns if not read.is_unmapped and bamFile.references[read.tid] in chrom_list['human']])+len([read for read in q2aligns if not read.is_unmapped and bamFile.references[read.tid] in chrom_list['human']])) > 0 and (len([read for read in q2aligns if read.is_unmapped]) + len([read for read in q1aligns if read.is_unmapped])) > 0:
    unknownReads += 1
    for b in q1aligns + q2aligns:
      seq = Counter(q1aligns[0].seq + q2aligns[0].seq)
      if 'N' in seq and seq['N'] >= 5:
        continue
      unknownFile.write(b)
  if (miscFile is not None and len([read for read in q2aligns if read.is_unmapped]) + len([read for read in q1aligns if read.is_unmapped])) == 0:
    for b in q1aligns + q2aligns:
      seq = Counter(q1aligns[0].seq + q2aligns[0].seq)
      if 'N' in seq and seq['N'] >= 5:
        continue
      miscFile.write(b)

transFile.close()
bamFile.close()
viralFile.close()
unknownFile.close()
if miscFile is not None:
  miscFile.close()
print totalReads, transReads, viralReads
