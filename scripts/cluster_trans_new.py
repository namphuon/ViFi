#!/usr/bin/env python

#cluster reads if consective reads have start position difference of less than 300bp regardless of strand                       

import pysam
import argparse
from time import clock
from collections import defaultdict, Counter
from sets import Set
import re

import hg19util as hg

parser = argparse.\
ArgumentParser(description="Cluster reads with one end mapping to viral to indicate insertion")
parser.add_argument('--data', dest='dataName',
                    help="Coordinatesorted BAM file aligned to combined hg19_viral reference with one end mapping to human genome and other to viral genome", metavar='FILE',
                    action='store', type=str, nargs=1)
parser.add_argument('--output', dest='outputName',
                    help="Output text file of putative insertion positions", metavar='FILE',
                    action='store', type=str, nargs=1)
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

args = parser.parse_args()
bamFile = pysam.Samfile(args.dataName[0], 'rb')
outFile = open(args.outputName[0], 'w')
rangeFile = open(args.outputName[0] + ".range", 'w')
hg19refs = Set(map(lambda x: 'chr' + str(x), range(1,23) + ['X', 'Y', 'M']) + map(str, range(1,23) +  ['X', 'Y']))
if args.chrom_list is not None:
  input = open(args.chrom_list[0], 'r')
  foo = [hg19refs.add(l) for l in input.next().strip().split(' ')]

MIN_SUPPORT = 3

def find_true_breakpoint_range(reads):
  pattern = re.compile('(\d+)([A-Z])')
  human_viral = re.compile('HHHHHX?X?X?VVVVV')
  viral_human = re.compile('VVVVVX?X?X?HHHHH')
  
  forward = Counter()
  reverse = Counter()
  for read in reads:    
    is_split = False
    if read.has_tag('SA'):
      tags = read.get_tag('SA').strip().split(';')
      for tag in tags:
        if len(tag) == 0:
          continue
        res = tag.split(',')
        if res[0].find('chr') == -1:
          virus = "".join([b1 * int(a1) for (a1,b1) in pattern.findall(res[3]) if b1 != 'D'])
          human = "".join([b1 * int(a1) for (a1,b1) in pattern.findall(read.cigarstring) if b1 != 'D'])
          combined = ""
          for x in range(len(virus)):
            if human[x] =='M' and virus[x] != 'M':
              combined+='H'
            elif human[x] !='M' and virus[x] == 'M':
              combined+='V'
            else:
              combined+='X'
          rev_combined = ''    
          virus = virus[::-1]        
          for x in range(len(virus)):
            if human[x] =='M' and virus[x] != 'M':
              rev_combined+='H'
            elif human[x] !='M' and virus[x] == 'M':
              rev_combined+='V'
            else:
              rev_combined+='X'
          if len(human_viral.findall(rev_combined)) != 0 or len(viral_human.findall(rev_combined)) != 0 or len(human_viral.findall(combined)) != 0 or len(viral_human.findall(combined)) != 0:
            is_split = True
    if read.is_reverse:
      reverse.update(read.positions)
      if is_split == True:
        reverse[read.positions[0]]+=10000      
    else:
      forward.update(read.positions)
      if is_split == True:
        forward[read.positions[-1]]+=10000      
    max_forward = forward.items()    
    max_forward.sort()
    max_reverse = reverse.items()
    max_reverse.sort()
    ranges = [-1, -1]
    splits = [-1, -1]
    if len(max_forward) >= 1:
      ranges[1] = max_forward[-1][0]
    if len(max_reverse) >= 1:
      ranges[0] = max_reverse[0][0]
    if ranges[0] == -1:
      ranges[0] = ranges[1]-300
    if ranges[1] == -1:
      ranges[1] = ranges[0]+300        
    if len(forward) > 0 and forward.most_common()[0][1] >= 10000:
      splits[0] = forward.most_common()[0][0]
    if len(reverse) > 0 and reverse.most_common()[0][1] >= 10000:
      splits[1] = reverse.most_common()[0][0]
  return((min(ranges),max(ranges), splits[0], splits[1]))
# find breakpoint that defines largest subset of alignments such that all forward alignments occur before all reverse alignments and at least one forward and on reverse alignment is included
# not used
def largest_clean_subset(clist):
        nr = len([a for a in clist if a.is_reverse])
        nf = len(clist) - nr
        maxpos = -1
        mr = 0
        mf = 0
        rc = nr
        fc = 0
        for a in clist:
                if not a.is_reverse:
                        fc += 1
                if fc >-1  and rc > -1 and (rc + fc > mr + mf):
                        maxpos = a.pos
                        mr = rc
                        mf = fc
                if a.is_reverse:
                        rc -= 1
        return (maxpos, mr, mf)

#takes a list of bam alignments and returns True if they form a good cluster
#at least 4 alignments
#all alignments such that forward reads appear before reverse reads
def clean_genomic_cluster (clist):
        fmax = -1
        rmin = -1
#        ls = largest_clean_subset(clist)
#        if ls[1] + ls[2] >MIN_SUPPORT:
        if len(clist) > MIN_SUPPORT and len(Set([a.qname for a in clist])) > MIN_SUPPORT:
                return True
        return False

#cluster reads if consective reads have start position difference of less than 300bp regardless of strand                       
clusterList = []
clist = []
caln = None
vlist = defaultdict(lambda: [], {})
vreads = defaultdict(lambda: Set([]), {})

for a in bamFile:
        vlist[(a.qname, a.is_read1)].append(a)
        if a.tid == -1 or bamFile.getrname(a.tid) not in hg19refs:
                if a.tid == -1:
                    continue
                vreads[bamFile.getrname(a.tid)].add((a.qname, a.is_read1))
                continue
        if caln is not None and (a.pos > caln.pos + 300 or caln.tid != a.tid) and clean_genomic_cluster(clist):
                clusterList.append(clist)
        if caln is None or a.pos > caln.pos + 300 or caln.tid != a.tid:
                clist = []
        caln = a
#        if hg.interval(a, bamfile=bamFile).num_unmasked() >= 35:
        if hg.interval(a, bamfile=bamFile).rep_content() <= 3 and a.mapq >= 10:
            clist.append(a)
if caln is not None and (a.pos > caln.pos + 300 or caln.tid != a.tid) and clean_genomic_cluster(clist):
        clusterList.append(clist)

clusterList.sort(lambda x, y: hg.interval(bamFile.getrname(x[0].tid), x[0].pos, x[-1].pos + x[-1].infer_query_length()) > hg.interval(bamFile.getrname(y[0].tid), y[0].pos, y[-1].pos + y[-1].infer_query_length()))

vsuper = {v: Set([v2 for v2 in vreads if v2 != v and len(vreads[v]) < len(vreads[v2]) and vreads[v].issubset(vreads[v2])]) for v in vreads}
vequaldict = {}
vequal = []
for v in vreads:
        inserted = False
        for vset in vequal:
                if vreads[v] == vreads[vset[0]]:
                        vset[1].add
                        vequaldict[v] = vset
                        inserted = True
                        break
        if not inserted:
                vset = (v, Set([v]))
                vequal.append(vset)
                vequaldict[v] = vset
rangeFile.write('Chr,Min,Max,Split1,Split2\n')
#outFile.write('#chr\tminpos\tmaxpos\t#reads\t#forward\t#reverse\t#viruses\tvirusnames\n')
outFile.write('#chr\tminpos\tmaxpos\t#reads\t#forward\t#reverse\n')
frbin = defaultdict(lambda: 0, {})
for x in range(21):
        for y in range(21):
                frbin[(x,y)] = 0

clusterSets = []
for c in clusterList:
        clusterSets.append(Set([a.qname for a in c]))

intersectionGraph = defaultdict(lambda: [], {})
cluster_duke35 = [hg.interval(bamFile.getrname(c[0].tid), c[0].pos, c[-1].pos + c[-1].infer_query_length()).rep_content()
                                    for c in clusterList]

for ci in range(len(clusterList)):
        c = clusterList[ci]
        (m1,m2,s1,s2) = find_true_breakpoint_range(c)    
        rangeFile.write("%s,%d,%d,%d,%d\n" % (bamFile.getrname(c[0].tid), m1, m2, s1, s2))
    
        vcount = defaultdict(lambda: [], {})
        for a in c:
                for a2 in vlist[(a.qname, not a.is_read1)]:
                        vcount[bamFile.getrname(a2.tid)].append(a2)
        vrep = Set([v for v in vcount.keys() if len(Set([a.qname for a in vcount[v]])) > MIN_SUPPORT
                    and len([a for a in vcount[v] if a.is_reverse]) > -1 and len([a for a in vcount[v] if not a.is_reverse]) > -1
                   ]
                  )
        cOvlList = []
        cInList = []
        cOutList = []
        maxOvl = 0
        maxOvlc = [None, -1, -1, 0]
        for c2i in range(len(clusterList)):
                if ci == c2i: continue
                if len(clusterSets[ci].intersection(clusterSets[c2i])) > maxOvl:
                        maxOvl = len(clusterSets[ci].intersection(clusterSets[c2i]))
                        c2 = clusterList[c2i]
                        maxOvlc = ([bamFile.getrname(c2[0].tid), c2[0].pos, c2[-1].pos+c2[-1].infer_query_length()])
                if clusterSets[ci].issubset(clusterSets[c2i]) and len(clusterSets[ci]) < len(clusterSets[c2i]):
                        c2 = clusterList[c2i]
                        cOutList.append([bamFile.getrname(c2[0].tid), c2[0].pos, c2[-1].pos+c2[-1].infer_query_length()])
                elif clusterSets[ci].issuperset(clusterSets[c2i]) and len(clusterSets[ci]) > len(clusterSets[c2i]):
                        c2 = clusterList[c2i]
                        cInList.append([bamFile.getrname(c2[0].tid), c2[0].pos, c2[-1].pos+c2[-1].infer_query_length()])
                elif len(clusterSets[ci].intersection(clusterSets[c2i])) > 0.9 * len(clusterSets[ci]) and len(clusterSets[ci].intersection(clusterSets[c2i])) > 0.9 * len(clusterSets[c2i]): 
                        c2 = clusterList[c2i]
                        cOvlList.append([bamFile.getrname(c2[0].tid), c2[0].pos, c2[-1].pos+c2[-1].infer_query_length()])
                elif len(clusterSets[ci].intersection(clusterSets[c2i])) > 0.9 * len(clusterSets[ci]): 
                        c2 = clusterList[c2i]
                        cOutList.append([bamFile.getrname(c2[0].tid), c2[0].pos, c2[-1].pos+c2[-1].infer_query_length()])
                elif len(clusterSets[ci].intersection(clusterSets[c2i])) > 0.9 * len(clusterSets[c2i]): 
                        c2 = clusterList[c2i]
                        cInList.append([bamFile.getrname(c2[0].tid), c2[0].pos, c2[-1].pos+c2[-1].infer_query_length()])
        ls = largest_clean_subset(c)
        cs = (len([a for a in c if not a.is_reverse]), len([a for a in c if a.is_reverse]))
        vplist = [(v, len(vcount[v]), min([a2.pos for a2 in vcount[v] if a2.is_reverse]+[100000000000]), max([-1]+[a2.pos+(a2.infer_query_length() if a2.infer_query_length() is not None else a2.qlen) for a2 in vcount[v] if not a2.is_reverse])) for v in vrep]
        vplist.sort(lambda x, y: y[1] - x[1])
        frbin[cs] += 1
        outFile.write("##==========================================================================================================================================================================================================================\n")
        #outFile.write('\t'.join(map(str, [bamFile.getrname(c[0].tid), c[0].pos, c[-1].pos + c[-1].infer_query_length(), len(Set([a.qname for a in c])), cs[0], cs[1], ls[1], ls[2], len([v for v in vcount if len(vcount[v]) > 3]), len(vrep), len(cOvlList), len(cInList), len(cOutList), 1/cluster_duke35[ci], maxOvl, maxOvlc] + vplist + cOvlList + cInList + cOutList)) + '\n')   
        outFile.write('\t'.join(map(str, [bamFile.getrname(c[0].tid), c[0].pos, c[-1].pos + c[-1].infer_query_length(), len(Set([a.qname for a in c])), cs[0], cs[1]])) + '\n')
        for a in c:
                outFile.write('##' + '\t'.join(map(str, [a.qname, bamFile.getrname(a.tid), a.pos, not a.is_reverse, a.is_read1])) + '\n')
        for v in vcount:
                for a2 in vcount[v]:
                        outFile.write('##' + '\t'.join(map(str, [a2.qname, bamFile.getrname(a2.tid), a2.pos, not a2.is_reverse, a2.is_read1])) + '\n')

outFile.close()
bamFile.close()
rangeFile.close()
print len(clusterList)

#for x in range(21):
#        for y in range(21):
#                cs = (x, y)
#                print x, y, frbin[cs]
#                print x, y + 1, frbin[cs]
#        print
#        for y in range(21):
#                cs = (x, y)
#                print x + 1, y, frbin[cs]
#                print x + 1, y + 1, frbin[cs]
#        print