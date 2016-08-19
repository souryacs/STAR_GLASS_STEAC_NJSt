import Header
from Header import *

#---------------------------------------------
"""
the function is used for checking equality of two floating point numbers
"""
def FlEq(a, b, eps=0.000001):
	#return (abs(math.log(a) - math.log(b)) <= eps)
	return (abs(a - b) <= eps)

#--------------------------------------------------
"""
this function returns the label of an internal or a leaf node 
in terms of newick representation
"""
def Node_Label(inp_node):
	return str(inp_node.as_newick_string(suppress_edge_lengths=True))

#-----------------------------------------------------
"""
this function reads the input tree list file
@parameters:
ROOTED_TREE and PRESERVE_UNDERSCORE are used for the dendropy based tree reading
INPUT_FILE_FORMAT: depending on its value, the trees are read from the input file 
according to NEWICK or NEXUS format
INPUT_FILENAME: the input file containing the treelist
"""
def Read_Gene_Data_Treelist(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME):
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
# this function defines couplet internode count with respect to the MRCA between two nodes
#--------------------------------------------------------
def DefineAccBranch(lca_node_level, node1, node2, GENE_TREE_TYPE):

	sum_of_branch_count = ((node1.level() - lca_node_level) + (node2.level() - lca_node_level)) - 1

	if (GENE_TREE_TYPE == MULTI_ALLELE):
		key1 = (AlleleToTaxon(node1.taxon.label), AlleleToTaxon(node2.taxon.label))
		key2 = (AlleleToTaxon(node2.taxon.label), AlleleToTaxon(node1.taxon.label))
		#print 'key1: ', key1, 'key2: ', key2
	else:
		key1 = (node1.taxon.label, node2.taxon.label)
		key2 = (node2.taxon.label, node1.taxon.label)

	""" 
	we insert the following information for this couplet
	1) increment the number of supporting tree
	2) add the couplet branch count value
	"""
	if key1 in TaxaPair_Reln_Dict:
		TaxaPair_Reln_Dict[key1]._IncrSupportTreeCount()
		TaxaPair_Reln_Dict[key1]._AddLevel(sum_of_branch_count)
	elif key2 in TaxaPair_Reln_Dict:
		TaxaPair_Reln_Dict[key2]._IncrSupportTreeCount()
		TaxaPair_Reln_Dict[key2]._AddLevel(sum_of_branch_count)
	else:
		TaxaPair_Reln_Dict.setdefault(key1, Reln_TaxaPair())
		TaxaPair_Reln_Dict[key1]._IncrSupportTreeCount()
		TaxaPair_Reln_Dict[key1]._AddLevel(sum_of_branch_count)

	return

#--------------------------------------------------------
# this function defines coalescence time with respect to the MRCA between two nodes
#--------------------------------------------------------
def DefineCoalTime(mrca_node, node1, node2, GENE_TREE_TYPE):
  
	coal_time1 = node1.distance_from_root() - mrca_node.distance_from_root()
	coal_time2 = node2.distance_from_root() - mrca_node.distance_from_root()
	sum_coal_time = coal_time1 + coal_time2
		
	if (GENE_TREE_TYPE == MULTI_ALLELE):
		key1 = (AlleleToTaxon(node1.taxon.label), AlleleToTaxon(node2.taxon.label))
		key2 = (AlleleToTaxon(node2.taxon.label), AlleleToTaxon(node1.taxon.label))
		#print 'key1: ', key1, 'key2: ', key2
	else:
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
def DefineCoalRank(lca_node_rank, node1, node2, GENE_TREE_TYPE):
	if (GENE_TREE_TYPE == MULTI_ALLELE):
		key1 = (AlleleToTaxon(node1.taxon.label), AlleleToTaxon(node2.taxon.label))
		key2 = (AlleleToTaxon(node2.taxon.label), AlleleToTaxon(node1.taxon.label))
		#print 'key1: ', key1, 'key2: ', key2
	else:
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
# this function derives couplet relations belonging to one tree
# that is provided as an input argument to this function
#--------------------------------------------------------
def DeriveCoupletRelations(Curr_tree, METHOD_USED, GENE_TREE_TYPE):
	# number of taxa in the current tree
	no_of_taxa = len(Curr_tree.infer_taxa().labels())
		
	# traverse the internal nodes of the tree in postorder fashion
	for curr_node in Curr_tree.postorder_internal_node_iter():
		# compute the coalescence rank associated with this node
		# this rank value will be used for all taxa pairs underlying this node
		curr_node_rank = no_of_taxa - curr_node.level()
		curr_node_level = curr_node.level()
					
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
					if (GENE_TREE_TYPE == MULTI_ALLELE):
						"""
						this check is required for multi allelle gene trees
						here leaf nodes may associate identical labels
						"""
						if (Same_Taxon(curr_node_child_leaf_nodes[i], curr_node_child_leaf_nodes[j])):
							#print 'same taxon (allelle) ', curr_node_child_leaf_nodes[i].taxon.label, ' and ', curr_node_child_leaf_nodes[j].taxon.label
							continue  
					"""
					otherwise we process this couplet containing different taxa information
					"""
					if (METHOD_USED == STAR):
						DefineCoalRank(curr_node_rank, curr_node_child_leaf_nodes[i], curr_node_child_leaf_nodes[j], GENE_TREE_TYPE)  
					elif (METHOD_USED == GLASS) or (METHOD_USED == STEAC):
						DefineCoalTime(curr_node, curr_node_child_leaf_nodes[i], curr_node_child_leaf_nodes[j], GENE_TREE_TYPE)
					else:
						DefineAccBranch(curr_node_level, curr_node_child_leaf_nodes[i], curr_node_child_leaf_nodes[j], GENE_TREE_TYPE)
			
		# one leaf node (direct descendant) and another leaf node (under one internal node)
		# will be related by ancestor / descendant relations
		if (len(curr_node_child_leaf_nodes) > 0) and (len(curr_node_child_internal_nodes) > 0):
			for p in curr_node_child_leaf_nodes:
				for q in curr_node_child_internal_nodes:
					for r in q.leaf_nodes():
						if (GENE_TREE_TYPE == MULTI_ALLELE):
							"""
							this check is required for multi allelle gene trees
							here leaf nodes may associate identical labels
							"""
							if (Same_Taxon(p, r)):
								#print 'same taxon (allelle) ', p.taxon.label, ' and ', r.taxon.label
								continue  
						"""
						otherwise we process this couplet containing different taxa information
						"""
						if (METHOD_USED == STAR):
							DefineCoalRank(curr_node_rank, p, r, GENE_TREE_TYPE)
						elif (METHOD_USED == GLASS) or (METHOD_USED == STEAC):
							DefineCoalTime(curr_node, p, r, GENE_TREE_TYPE)
						else:
							DefineAccBranch(curr_node_level, p, r, GENE_TREE_TYPE)
		
		# finally a pair of leaf nodes which are descendant of internal nodes will be related by NO_EDGE relation
		if (len(curr_node_child_internal_nodes) > 1):
			for i in range(len(curr_node_child_internal_nodes) - 1):
				for j in range(i+1, len(curr_node_child_internal_nodes)):
					for p in curr_node_child_internal_nodes[i].leaf_nodes():
						for q in curr_node_child_internal_nodes[j].leaf_nodes():
							if (GENE_TREE_TYPE == MULTI_ALLELE):
								"""
								this check is required for multi allelle gene trees
								here leaf nodes may associate identical labels
								"""
								if (Same_Taxon(p, q)):
									#print 'same taxon (allelle) ', p.taxon.label, ' and ', q.taxon.label
									continue  
							"""
							otherwise we process this couplet containing different taxa information
							"""      
							if (METHOD_USED == STAR):
								DefineCoalRank(curr_node_rank, p, q, GENE_TREE_TYPE)
							elif (METHOD_USED == GLASS) or (METHOD_USED == STEAC):
								DefineCoalTime(curr_node, p, q, GENE_TREE_TYPE)
							else:
								DefineAccBranch(curr_node_level, p, q, GENE_TREE_TYPE)

	return

#-----------------------------------------------------
# this function converts an allele to their main taxon label
#-----------------------------------------------------
def AlleleToTaxon(node_label):
	if node_label in Allele_Dict:
		return Allele_Dict[node_label]
	return node_label
  
#-----------------------------------------------------
# this function checks if the input taxa pair belong from the same taxon (they are allelles)
# i.e. their labels are either equivalent or from the same taxon 
#-----------------------------------------------------
def Same_Taxon(node1, node2):
	node1_label = node1.taxon.label
	node2_label = node2.taxon.label
	if (AlleleToTaxon(node1_label) == AlleleToTaxon(node2_label)):
		return True
	return False
  
##-----------------------------------------------------
## this function finds the MRCA of this two input taxa labels
##-----------------------------------------------------
#def Find_MRCA(Inp_Tree, spec_list):
  #node1 = Inp_Tree.find_node_with_taxon_label(spec_list[0])
  #pn = node1.parent_node
  #while (pn is not None):
    #leaf_labels = []
    #for n in pn.leaf_nodes():
      #leaf_labels.append(n.taxon.label)
    #if set(spec_list).issubset(set(leaf_labels)):
      #return pn
    #pn = pn.parent_node
      
  #return None
  
#-----------------------------------------------------
""" 
this function modifies input treelist file containing multiple alleles 
and generates allele information 
"""
#-----------------------------------------------------  
def ModifyMultiAlleleTreelist(inp_filename, modified_filename, Output_Text_File):

	"""
	this is the count of number of input phylogenetic trees
	"""
	tree_count = 0

	"""
	first we read all the lines of the original treelist file
	containing allelle information
	"""
	fp1 = open(inp_filename, 'r')
	lines = fp1.readlines()
	fp1.close()
	fp2 = open(modified_filename, 'w')
	for lc in range(len(lines)):
		curr_line = lines[lc]
		"""
		from the current line, delete all leading trailing spaces and newline characters
		"""
		curr_line.strip('')
		curr_line = re.sub('[\n]', '', curr_line)
		if (curr_line[len(curr_line) - 1] == ';'):
			""" 
			this line represents a tree string
			it may contain multi allele tree
			"""
			#print 'original line: ', curr_line
			
			# increment the tree counter
			tree_count = tree_count + 1
			
			#print '**** Tree count: ', tree_count
			
			# create a key of the global tree allele dictionary
			# initially an empty list is created
			Tree_Allele_Dict.setdefault((tree_count - 1), [])
			
			# remove rooted, unrooted information
			curr_line = re.sub('[\[&U\]]', '', curr_line)
			curr_line = re.sub('[\[&u\]]', '', curr_line)
			curr_line = re.sub('[\[&R\]]', '', curr_line)
			curr_line = re.sub('[\[&r\]]', '', curr_line)
			# remove branch length information
			curr_line = re.sub('[:][0-9]*.[0-9]*', '', curr_line)
			# remove brackets and semicolon
			curr_line = re.sub('[();]', '', curr_line)
			#print 'after re sub -  line: ', curr_line
			"""
			now we insert individual taxon information and their frequencies (>1 means allelle)
			in two different custom defined lists
			"""
			taxon_list = []
			freq_list = []
			for t in curr_line.split(','):
				if (len(taxon_list) == 0):
					taxon_list.append(t)
					freq_list.append(1)
				elif t not in taxon_list:
					taxon_list.append(t)
					freq_list.append(1)
				else:
					idx = taxon_list.index(t)
					freq_list[idx] = freq_list[idx] + 1

			#print 'taxon_list: ', taxon_list
			#print 'freq_list: ', freq_list
			
			"""
			we scan this derived taxon list and insert it in the COMPLETE_INPUT_TAXA_LIST
			"""
			for t in taxon_list:
				if t not in COMPLETE_INPUT_TAXA_LIST:
					COMPLETE_INPUT_TAXA_LIST.append(t)
			
			""" 
			now we have to modify original line with this taxon information
			we scan through the taxon and frequency lists and check multiple frequency taxon
			"""
			# this is the original treelist content
			orig_line = lines[lc]
			#print 'orig_line: ', orig_line
			for freq_idx in range(len(freq_list)):
				f = freq_list[freq_idx]
				if (f > 1):
					""" 
					current taxon has multiple allelles within this tree string
					"""
					t = taxon_list[freq_idx]
					#print 'taxon with multiple occurrence: ', t
					
					"""
					insert the taxon information and its frequency count in the tree allele dictionary
					the sublist has following three elements:
					1) taxon name whose multi allele information is present
					2) number of allele instance
					3) index of taxon which has minimum level information - initially it is set as 0 - it will be computed later
					"""
					subl = [t, f, 0]
					key1 = (tree_count - 1)
					Tree_Allele_Dict[key1].append(subl)
					
					"""
					when searching for this taxon name within the original tree encoded string,
					this taxon name may be errorneously embedded within a bigger taxon name
					for example, taxon_1 may be embedded within taxon_11 
					so, we check for the character : or , after the taxon name
					for unweighted trees, a taxon name may b followed by a closing first bracket ) as well
					so we use "\)" characters where \ is the escape character
					"""
					"""
					we note all the start and end positions of taxon name occurrence
					"""
					t1 = t + ':'
					t2 = t + ','
					t3 = t + ')'
					
					#------------------------------------------------------
					# comment - sourya
					#all_start_pos_list = [[m.start(), m.end()] for m in re.finditer(t1, orig_line)]
					#all_start_pos_list.extend([[m.start(), m.end()] for m in re.finditer(t2, orig_line)])
					#all_start_pos_list.extend([[m.start(), m.end()] for m in re.finditer(t3, orig_line)])
					#print 'all_start_pos_list: ', all_start_pos_list
					#for idx in range(len(all_start_pos_list) - 1, 0, -1):
						#"""
						#we scan from the last of the tree string to the first
						#and replace the taxon occurrence with our custom defined allelle names
						#to make a string for parsing in dendropy package
						#"""
						## this is the allelle name (custom defined) corresponding to the current taxon
						#new_taxon_name = t + '_' + str(idx)
						#"""
						#insert this allelle information in a custom defined allelle dictionary
						#"""
						#if new_taxon_name not in Allele_Dict:
							#Allele_Dict.setdefault(new_taxon_name, t)
						#print 'new_taxon_name: ', new_taxon_name
						
						#"""
						#replace the taxon name by this allelle name in the original string representation of tree
						#"""
						#orig_line = orig_line[:all_start_pos_list[idx][0]] + new_taxon_name + orig_line[(all_start_pos_list[idx][1] - 1):]
					# end comment - sourya
					#------------------------------------------------------
					# add - sourya
					
					pattern_list = [t1, t2, t3]
					allele_count = f
					#print 'allele count: ', allele_count
					for tp in pattern_list:
						#print 'tp: ', tp
						while (allele_count > 1):
							"""
							we scan from the last of the tree string to the first
							and replace the taxon occurrence with our custom defined allelle names
							to make a string for parsing in dendropy package
							"""
							idx = orig_line.rfind(tp)
							if (idx == -1):
								break
							#print 'pattern : ', tp, ' is found at idx: ', idx
							allele_count = allele_count - 1
							if (allele_count > 0):
								# t is the original taxon name
								# this is the allelle name (custom defined) corresponding to the current taxon
								new_taxon_name = t + '_' + str(allele_count)
								"""
								insert this allelle information in a custom defined allelle dictionary
								"""
								if new_taxon_name not in Allele_Dict:
									Allele_Dict.setdefault(new_taxon_name, t)
								#print 'new_taxon_name: ', new_taxon_name
								"""
								replace the taxon name by this allelle name in the original string representation of tree
								"""
								orig_line = orig_line[:idx] + new_taxon_name + orig_line[(idx + len(t)):]
								#print 'orig line - middle phase: ', orig_line
					
					# end add - sourya
					#------------------------------------------------------
						
			#print 'orig_line after taxon replace: ', orig_line
			"""
			write the modified tree in the output file
			"""
			fp2.write(str(orig_line))
		else:  
			# this string does not correspond to a tree string
			# copy the string as it is
			fp2.write(str(curr_line))    
			# if this line (with respect to input file) is not the last line
			# then append a newline character
			if (lc < (len(lines) - 1)):
				fp2.write('\n')

	# close the output file
	fp2.close()

	"""
	write the allelle dictionary in the output text file
	"""
	if (DEBUG_LEVEL > 2):
		fp3 = open(Output_Text_File, 'a')
		fp3.write('\n\n This input file contains gene trees with multi allelle information')
		fp3.write('\n\n Below we write the allelles along with the original taxon name information')
		fp3.write('\n\n COMPLETE_INPUT_TAXA_LIST: ' + str(COMPLETE_INPUT_TAXA_LIST))
		fp3.write('\n\n Allele_Dict: ' + str(Allele_Dict))
		for tr_key in Tree_Allele_Dict:
			fp3.write('\n Tree allele dict key: ' + str(tr_key) + ' contents: ' + str(Tree_Allele_Dict[tr_key]))
		fp3.close()

	return
  