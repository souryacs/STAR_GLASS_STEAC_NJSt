import Header
from Header import *

#--------------------------------------------------
# this function returns the label of an internal or a leaf node 
# in terms of newick representation
def Node_Label(inp_node):
  return str(inp_node.as_newick_string(suppress_edge_lengths=True))

##-----------------------------------------------------
''' this function reads the input tree collection file
the file contains a collection of input candidate source trees
each such tree is composed of a large no of taxa (placed at the leaves of the tree) '''
def Read_Gene_Data_Treelist(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME):
  ''' depending on the value of INPUT_FILE_FORMAT
  the data is read from the file according to NEWICK or NEXUS format '''
  if (INPUT_FILE_FORMAT == 1):
    Gene_TreeList = dendropy.TreeList.get_from_path(INPUT_FILENAME, \
						      schema="newick", \
						      preserve_underscores=PRESERVE_UNDERSCORE, \
						      default_as_rooted=ROOTED_TREE)
  else:
    Gene_TreeList = dendropy.TreeList.get_from_path(INPUT_FILENAME, \
						      schema="nexus", \
						      preserve_underscores=PRESERVE_UNDERSCORE, \
						      default_as_rooted=ROOTED_TREE)
  
  return Gene_TreeList

#--------------------------------------------------------
# this function defines coalescence time with respect to the MRCA between two nodes
#--------------------------------------------------------
def DefineCoalTime(mrca_node, node1, node2):
  
  coal_time1 = node1.distance_from_root() - mrca_node.distance_from_root()
  coal_time2 = node2.distance_from_root() - mrca_node.distance_from_root()
  sum_coal_time = coal_time1 + coal_time2
    
  key1 = (node1.taxon.label, node2.taxon.label)
  key2 = (node2.taxon.label, node1.taxon.label)
  """ 
  we insert the following information for this couplet
  1) increment the number of supporting tree
  2) add the sum of coalescence time value with respect to the MRCA node between this couplet 
  """
  if key1 in TaxaPair_Reln_Dict:
    TaxaPair_Reln_Dict[key1]._IncrSupportTreeCount()
    TaxaPair_Reln_Dict[key1]._AddCoalescenceInfo(sum_coal_time)
  elif key2 in TaxaPair_Reln_Dict:
    TaxaPair_Reln_Dict[key2]._IncrSupportTreeCount()
    TaxaPair_Reln_Dict[key2]._AddCoalescenceInfo(sum_coal_time)
  else:
    TaxaPair_Reln_Dict.setdefault(key1, Reln_TaxaPair())
    TaxaPair_Reln_Dict[key1]._IncrSupportTreeCount()
    TaxaPair_Reln_Dict[key1]._AddCoalescenceInfo(sum_coal_time)
  
  return
      
#--------------------------------------------------------
# this function defines coalescence rank of the MRCA between two nodes
#--------------------------------------------------------
def DefineCoalRank(lca_node_rank, node1, node2):
  key1 = (node1.taxon.label, node2.taxon.label)
  key2 = (node2.taxon.label, node1.taxon.label)
  """ 
  we insert the following information for this couplet
  1) increment the number of supporting tree
  2) add the LCA rank value
  """
  if key1 in TaxaPair_Reln_Dict:
    TaxaPair_Reln_Dict[key1]._IncrSupportTreeCount()
    TaxaPair_Reln_Dict[key1]._AddRankInfo(lca_node_rank)
  elif key2 in TaxaPair_Reln_Dict:
    TaxaPair_Reln_Dict[key2]._IncrSupportTreeCount()
    TaxaPair_Reln_Dict[key2]._AddRankInfo(lca_node_rank)
  else:
    TaxaPair_Reln_Dict.setdefault(key1, Reln_TaxaPair())
    TaxaPair_Reln_Dict[key1]._IncrSupportTreeCount()
    TaxaPair_Reln_Dict[key1]._AddRankInfo(lca_node_rank)

  return
      
#--------------------------------------------------------
# this function derives coupket relations belonging to one tree
# that is provided as an input argument to this function
#--------------------------------------------------------
def DeriveCoupletRelations(Curr_tree, METHOD_USED):
  # number of taxa in the current tree
  no_of_taxa = len(Curr_tree.infer_taxa().labels())
    
  # traverse the internal nodes of the tree in postorder fashion
  for curr_node in Curr_tree.postorder_internal_node_iter():
    # compute the coalescence rank associated with this node
    # this rank value will be used for all taxa pairs underlying this node
    curr_node_rank = no_of_taxa - curr_node.level()
          
    # list the leaf and internal children of the current node
    curr_node_child_leaf_nodes = []
    curr_node_child_internal_nodes = []
    for x in curr_node.child_nodes():
      if (x.is_leaf() == True):
	curr_node_child_leaf_nodes.append(x)
      else:
	curr_node_child_internal_nodes.append(x)
    
    # pair of leaf nodes will be related by sibling relations
    if (len(curr_node_child_leaf_nodes) > 1):
      for i in range(len(curr_node_child_leaf_nodes) - 1):
	for j in range(i+1, len(curr_node_child_leaf_nodes)):
	  if (METHOD_USED == STAR):
	    DefineCoalRank(curr_node_rank, curr_node_child_leaf_nodes[i], curr_node_child_leaf_nodes[j])  
	  else:
	    DefineCoalTime(curr_node, curr_node_child_leaf_nodes[i], curr_node_child_leaf_nodes[j])
	    
    # one leaf node (direct descendant) and another leaf node (under one internal node)
    # will be related by ancestor / descendant relations
    if (len(curr_node_child_leaf_nodes) > 0) and (len(curr_node_child_internal_nodes) > 0):
      for p in curr_node_child_leaf_nodes:
	for q in curr_node_child_internal_nodes:
	  for r in q.leaf_nodes():
	    if (METHOD_USED == STAR):
	      DefineCoalRank(curr_node_rank, p, r)
	    else:
	      DefineCoalTime(curr_node, p, r)
    
    # finally a pair of leaf nodes which are descendant of internal nodes will be related by NO_EDGE relation
    if (len(curr_node_child_internal_nodes) > 1):
      for i in range(len(curr_node_child_internal_nodes) - 1):
	for j in range(i+1, len(curr_node_child_internal_nodes)):
	  for p in curr_node_child_internal_nodes[i].leaf_nodes():
	    for q in curr_node_child_internal_nodes[j].leaf_nodes():
	      if (METHOD_USED == STAR):
		DefineCoalRank(curr_node_rank, p, q)
	      else:
		DefineCoalTime(curr_node, p, q)
  
  return
  
#-----------------------------------------------------
# this function finds the MRCA of this two input taxa labels
#-----------------------------------------------------
def Find_MRCA(Inp_Tree, spec_list):
  node1 = Inp_Tree.find_node_with_taxon_label(spec_list[0])
  pn = node1.parent_node
  while (pn is not None):
    leaf_labels = []
    for n in pn.leaf_nodes():
      leaf_labels.append(n.taxon.label)
    if set(spec_list).issubset(set(leaf_labels)):
      return pn
    pn = pn.parent_node
      
  return None
  
  
