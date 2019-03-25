#!/usr/bin/env python

#Builds a collection of HMMs from an existing alignment and tree

import argparse, pysam, copy, dendropy, os
from dendropy import Tree

def decompose_tree(tree, max_size = 10, tree_map = {}, decomposition = 'hierarchical' ):
  """
  This function decomposes the tree until all subtrees are smaller than 
  the max size, but does not decompose below min size.  
  Returns a map containing the subtrees, in an ordered fashion.  
  """          
  #Don't deroot if doing clade-based decomposition
  tree.deroot()  
  nodes = len(tree.leaf_nodes())*1.0  
  if (decomposition == 'standard' and nodes <= max_size) or decomposition == 'hierarchical':
    tree_map[len(tree_map)] = copy.deepcopy(tree)
  if nodes > max_size:
    splits = [(n, len(n.leaf_nodes())/nodes) for n in tree.nodes()]
    splits = sorted(splits, key = lambda s:abs(0.50-s[1]))
    node_filter_right = lambda nd: nd in splits[0][0].leaf_nodes() and nd.is_leaf()
    node_filter_left = lambda nd: nd not in splits[0][0].leaf_nodes() and nd.is_leaf()
    right = tree.extract_tree(node_filter_fn=node_filter_right)
    left = tree.extract_tree(node_filter_fn=node_filter_left)
    decompose_tree(right, max_size = max_size, tree_map=tree_map, decomposition=decomposition)
    decompose_tree(left, max_size = max_size, tree_map=tree_map, decomposition=decomposition)
  return tree_map

def build_hmms(tree_map, alignment_file, output_dir, prefix, keep_alignment):
  """
  This function takes the decomposed trees, builds a FASTA alignment from the 
  induced alignment from the leaves in the subtrees, and builds an HMM from
  the induced alignment
  """          
  alignment = dendropy.DnaCharacterMatrix.get(path=alignment_file, schema='fasta')  
  if not os.path.exists(output_dir):
    os.makedirs(output_dir)
  for tree in tree_map:    
    out_file_prefix = "%s/%s.%d" % (output_dir,prefix,tree)
    align_file = "%s.aln" % (out_file_prefix)
    names = [n.taxon.label for n in tree_map[tree].leaf_nodes()]
    pairs = [(k,v) for (k,v) in alignment.items() if k.label in names]
    output = open(align_file, 'w')
    for (k,v) in pairs:
      output.write('>%s\n%s\n' % (k.label, str(v)))
    output.close()
    os.system('hmmbuild %s.hmmbuild %s.aln' % (out_file_prefix, out_file_prefix))
    if not keep_alignment:
      os.remove('%s.aln' % out_file_prefix)

parser = argparse.ArgumentParser(description="Buids collection of HMMs from an alignment and tree")
parser.add_argument('--tree_file', dest='tree_file',
                    help="A NEWICK tree file containing the phylogenetic tree of the viral family of interest", metavar='TREE',
                    action='store', type=str)
parser.add_argument('--alignment_file', dest='alignment_file',
                    help="A FASTA file containing the alignment of the viral family of interest", metavar='ALIGNMENT',
                    action='store', type=str)                  
parser.add_argument('--output_dir', dest='output_dir',
                    help="output directory of the resulting HMMs", metavar='FILE',
                    action='store', type=str)
parser.add_argument('--prefix', dest='prefix',
                    help="prefix used to name HMM files (default=viral)", metavar='PREFIX', default = 'viral',
                    action='store', type=str)
parser.add_argument('--max_size_fraction', dest='max_size',
                    help="Maximum fraction of total size to decompose (default=.10)", metavar='MAX_SIZE', default = 0.10,
                    action='store', type=int)
parser.add_argument('--keep_alignment', dest='keep_alignment',
                    help="Keep temporary ", default = False,
                    action='store', type=bool)
parser.add_argument('--decomposition', dest='decomposition',
                    help="Decomposition strategy:standard, hierarchical (default=hierarchical)", default = 'hierarchical',
                    action='store', type=str)

arg = parser.parse_args()
tree = Tree.get(file=open(arg.tree_file, 'r'), schema="newick", preserve_underscores=True)
tree_map = {}
print "Decomposing Tree"
max_size = max(10, int(arg.max_size*len(tree.leaf_nodes())))
decompose_tree(tree, max_size, tree_map = tree_map, decomposition=arg.decomposition)
print "Building HMMs"
build_hmms(tree_map, arg.alignment_file, arg.output_dir, arg.prefix, arg.keep_alignment)