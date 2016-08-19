import Header
from Header import *
import UtilFunc
from UtilFunc import *

#---------------------------------------------
"""
this function refines the initial species tree (in terms of a star network) to 
find the true species tree
it does using agglomerative clustering (NJ principle)
the distance metric employed for NJ algorithm can vary depending on experimentation 
"""
def Form_Species_Tree_NJ_Cluster(Star_Tree_Initial, COMPLETE_INPUT_TAXA_LIST, METHOD_USED, Output_Text_File):
	"""
	initially we have N of clusters for N taxa, where individual clusters are isolated
	agglomerating technique introduces a bipartition (speciation) which contains two taxa as its children
	"""
	no_of_taxa_clust = len(COMPLETE_INPUT_TAXA_LIST)
	
	"""
	initialize the taxa clusters
	Note: the copy operation is done by individual elements
	"""
	clust_species_list = []
	for i in range(len(COMPLETE_INPUT_TAXA_LIST)):
		subl = []
		subl.append(COMPLETE_INPUT_TAXA_LIST[i])
		clust_species_list.append(subl)

	if (DEBUG_LEVEL >= 2):
		fp = open(Output_Text_File, 'a')
		fp.write('\n COMPLETE_INPUT_TAXA_LIST ' + str(COMPLETE_INPUT_TAXA_LIST))
		fp.write('\n Initial formed clust_species_list ' + str(clust_species_list))
		fp.close()        

	"""
	allocate a 2D square matrix of to store the pairwise distance measure for individual pairs of taxa clusters
	"""
	Dist_Mat_clust_pair_NJ = numpy.zeros((no_of_taxa_clust, no_of_taxa_clust), dtype=numpy.float)

	"""
	according to the species tree construction method, fill the "Dist_Mat_clust_pair_NJ" 
	with the appropriate distance measure
	"""
	for l in TaxaPair_Reln_Dict:
		spec1 = l[0]
		spec2 = l[1]
		spec1_idx = COMPLETE_INPUT_TAXA_LIST.index(spec1)
		spec2_idx = COMPLETE_INPUT_TAXA_LIST.index(spec2)
		if (METHOD_USED == STAR):
			# average coalescence rank
			Dist_Mat_clust_pair_NJ[spec1_idx][spec2_idx] = Dist_Mat_clust_pair_NJ[spec2_idx][spec1_idx] = TaxaPair_Reln_Dict[l]._GetAvgRank()
		elif (METHOD_USED == GLASS):
			# minimum coalescence time
			Dist_Mat_clust_pair_NJ[spec1_idx][spec2_idx] = Dist_Mat_clust_pair_NJ[spec2_idx][spec1_idx] = TaxaPair_Reln_Dict[l]._GetMinCoalescenceTime()
		elif (METHOD_USED == STEAC):
			# average coalescence time
			Dist_Mat_clust_pair_NJ[spec1_idx][spec2_idx] = Dist_Mat_clust_pair_NJ[spec2_idx][spec1_idx] = TaxaPair_Reln_Dict[l]._GetAvgCoalescenceTime()
		else:
			# average internode count
			Dist_Mat_clust_pair_NJ[spec1_idx][spec2_idx] = Dist_Mat_clust_pair_NJ[spec2_idx][spec1_idx] = TaxaPair_Reln_Dict[l]._GetAvgSumLevel()

	"""
	For NJ based agglomeration, 
	allocate one new square matrix which will contain the relative distances between individual cluster pairs
	"""
	Norm_DistMat_ClustPair_NJ = numpy.zeros((no_of_taxa_clust, no_of_taxa_clust), dtype=numpy.float)

	#-------------------------------------------------------
	"""
	loop to execute the agglomerative clustering
	"""
	while(no_of_taxa_clust > 2): 
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n iteration start --- number of clusters: ' + str(no_of_taxa_clust))
			fp.write('\n clust_species_list : ' + str(clust_species_list))
			fp.write('\n printing contents of Dist_Mat_clust_pair_NJ ---- ')
			for i in range(no_of_taxa_clust):
				fp.write('\n ' + str(i) + '--' + str(clust_species_list[i]) + '--->>')
				for j in range(i+1):
					fp.write(' ' + str(Dist_Mat_clust_pair_NJ[i][j]))     
			fp.close()
		
		"""
		for individual cluster Cx, it contains Dist(Cx, :) - sum of distance values (according to the distance metric employed) 
		for all the cluster pairs (Cx, Cy), where Cy is any other cluster
		"""
		sum_Dist_from_one_Clust_List = []
		for i in range(no_of_taxa_clust):
			t = 0
			for j in range(no_of_taxa_clust):
				t = t + Dist_Mat_clust_pair_NJ[i][j]
			sum_Dist_from_one_Clust_List.append(t)
			
		if (DEBUG_LEVEL > 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n content of sum_Dist_from_one_Clust_List : ' + str(sum_Dist_from_one_Clust_List))
			fp.close()
			
		"""
		fill the relative distance matrix "Norm_DistMat_ClustPair_NJ"
		individual elements: Hij = Xij - 1/(N-2)(Xi + Xj)
		where Xi = sum_Dist_from_one_Clust_List[i] - Dist_Mat_clust_pair_NJ[i][j]
		"""
		for i in range(no_of_taxa_clust - 1):
			for j in range(i+1, no_of_taxa_clust):
				"""
				here ri , rj are the sum of all distances
				"""
				ri = sum_Dist_from_one_Clust_List[i] / (no_of_taxa_clust - 2)
				rj = sum_Dist_from_one_Clust_List[j] / (no_of_taxa_clust - 2)
				"""
				relative distance matrix entries
				"""
				Norm_DistMat_ClustPair_NJ[i][j] = (Dist_Mat_clust_pair_NJ[i][j] - ri - rj)
				Norm_DistMat_ClustPair_NJ[j][i] = Norm_DistMat_ClustPair_NJ[i][j]

		"""
		printing the relative Distance matrix entries
		"""
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n printing contents of Norm_DistMat_ClustPair_NJ ---- ')
			for i in range(no_of_taxa_clust):
				fp.write('\n ' + str(i) + '--' + str(clust_species_list[i]) + '--->>')
				for j in range(i+1):
					fp.write(' ' + str(Norm_DistMat_ClustPair_NJ[i][j]))   
			fp.close()
		
		"""
		find the cluster pair having the minimum relative distance 
		(with respect to the matrix Norm_DistMat_ClustPair_NJ)
		"""
		min_val = Norm_DistMat_ClustPair_NJ[0][1]
		min_idx_i = 0
		min_idx_j = 1
		for i in range(no_of_taxa_clust - 1):
			for j in range(i+1, no_of_taxa_clust):
				if (i == j):
					continue
				if (Norm_DistMat_ClustPair_NJ[i][j] < min_val):
					min_val = Norm_DistMat_ClustPair_NJ[i][j]
					min_idx_i = i
					min_idx_j = j
				# add - sourya - 17/06/2016
				elif (FlEq(Norm_DistMat_ClustPair_NJ[i][j], min_val) == True):
					if (Dist_Mat_clust_pair_NJ[i][j] < Dist_Mat_clust_pair_NJ[min_idx_i][min_idx_j]):
						"""
						here the NJ based distance matrix values are same
						so, we check the original distance matrix entries
						and prefer the element with the lowest value
						"""
						min_idx_i = i
						min_idx_j = j
				# end add - sourya
				# comment - sourya
				# this heuristic is not a part of original STAR / STEAC / NJst algorithm
				#elif (Norm_DistMat_ClustPair_NJ[i][j] == min_val):
					## here we prioritize the cluster pair having minimum number of species
					#if (len(clust_species_list[i]) + len(clust_species_list[j])) < (len(clust_species_list[min_idx_i]) + len(clust_species_list[min_idx_j])):
						#min_idx_i = i
						#min_idx_j = j
				# end comment - sourya
		
		"""
		note down the taxa list in these two indices (min_idx_i and min_idx_j) of the clust_species_list
		"""
		taxa_list = []
		for x in clust_species_list[min_idx_i]:
			taxa_list.append(x)
		for x in clust_species_list[min_idx_j]:
			taxa_list.append(x)

		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n min_idx_i ' + str(min_idx_i) + ' min_idx_j : ' + str(min_idx_j) + 'min val : ' + str(min_val))
			fp.write('\n min_idx_i species list ' + str(clust_species_list[min_idx_i]))
			fp.write('\n min_idx_j species list ' + str(clust_species_list[min_idx_j]))
			fp.write('\n complete taxa list (union) ' + str(taxa_list))
			fp.close()

		#---------------------------------------------------------      
		"""
		*****   for individual clusters, we check if the cluster contains one or more species  ******
		"""
		if (len(clust_species_list[min_idx_i]) > 1) and (len(clust_species_list[min_idx_j]) > 1):
			"""
			case 1 - both the clusters have > 1 species
			and the clusters are represented by an internal node which is the MRCA of the constituent species set
			"""
			first_cluster_mrca_node = Star_Tree_Initial.mrca(taxon_labels=clust_species_list[min_idx_i])
			second_cluster_mrca_node = Star_Tree_Initial.mrca(taxon_labels=clust_species_list[min_idx_j])
			all_taxa_mrca_node = Star_Tree_Initial.mrca(taxon_labels=taxa_list)
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n label of first_cluster_mrca_node: ' + str(Node_Label(first_cluster_mrca_node)))      
				fp.write('\n label of second_cluster_mrca_node: ' + str(Node_Label(second_cluster_mrca_node)))
				fp.write('\n label of all_taxa_mrca_node: ' + str(Node_Label(all_taxa_mrca_node)))
				fp.close()

			# create new internal node 
			newnode = dendropy.Node()
			# its parent node will be the previous MRCA node of all the taxa in two clusters
			all_taxa_mrca_node.add_child(newnode)
			newnode.parent_node = all_taxa_mrca_node
			all_taxa_mrca_node.remove_child(first_cluster_mrca_node)
			first_cluster_mrca_node.parent_node = None
			all_taxa_mrca_node.remove_child(second_cluster_mrca_node)
			second_cluster_mrca_node.parent_node = None
			# add these individual clusters' MRCA node as its children
			newnode.add_child(first_cluster_mrca_node)
			first_cluster_mrca_node.parent_node = newnode
			newnode.add_child(second_cluster_mrca_node)
			second_cluster_mrca_node.parent_node = newnode
			# update splits of the resulting tree
			Star_Tree_Initial.update_splits(delete_outdegree_one=False)
		
		elif (len(clust_species_list[min_idx_i]) == 1) and (len(clust_species_list[min_idx_j]) > 1):
			"""
			case 2: second cluster has at least 2 species, while the first cluster is a leaf
			"""
			first_cluster_leaf_node = Star_Tree_Initial.find_node_with_taxon_label(clust_species_list[min_idx_i][0])
			second_cluster_mrca_node = Star_Tree_Initial.mrca(taxon_labels=clust_species_list[min_idx_j])
			all_taxa_mrca_node = Star_Tree_Initial.mrca(taxon_labels=taxa_list)
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n first cluster is a leaf - its label: ' + str(Node_Label(first_cluster_leaf_node)))      
				fp.write('\n label of second_cluster_mrca_node: ' + str(Node_Label(second_cluster_mrca_node)))
				fp.write('\n label of all_taxa_mrca_node: ' + str(Node_Label(all_taxa_mrca_node)))
				fp.close()
			
			# create new internal node 
			newnode = dendropy.Node()
			# its parent node will be the previous MRCA node of all the taxa in two clusters
			all_taxa_mrca_node.add_child(newnode)
			newnode.parent_node = all_taxa_mrca_node
			all_taxa_mrca_node.remove_child(first_cluster_leaf_node)
			first_cluster_leaf_node.parent_node = None
			all_taxa_mrca_node.remove_child(second_cluster_mrca_node)
			second_cluster_mrca_node.parent_node = None
			# add these individual clusters' MRCA node as its children
			newnode.add_child(first_cluster_leaf_node)
			first_cluster_leaf_node.parent_node = newnode
			newnode.add_child(second_cluster_mrca_node)
			second_cluster_mrca_node.parent_node = newnode
			# update splits of the resulting tree
			Star_Tree_Initial.update_splits(delete_outdegree_one=False)
			
		elif (len(clust_species_list[min_idx_i]) > 1) and (len(clust_species_list[min_idx_j]) == 1):
			"""
			case 3: first cluster has at least 2 species, while the second cluster is a leaf
			"""
			first_cluster_mrca_node = Star_Tree_Initial.mrca(taxon_labels=clust_species_list[min_idx_i])
			second_cluster_leaf_node = Star_Tree_Initial.find_node_with_taxon_label(clust_species_list[min_idx_j][0])
			all_taxa_mrca_node = Star_Tree_Initial.mrca(taxon_labels=taxa_list)
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n label of first_cluster_mrca_node: ' + str(Node_Label(first_cluster_mrca_node)))      
				fp.write('\n label of second_cluster_mrca_node: ' + str(Node_Label(second_cluster_leaf_node)))
				fp.write('\n label of all_taxa_mrca_node: ' + str(Node_Label(all_taxa_mrca_node)))
				fp.close()
			
			# create new internal node 
			newnode = dendropy.Node()
			# its parent node will be the previous MRCA node of all the taxa in two clusters
			all_taxa_mrca_node.add_child(newnode)
			newnode.parent_node = all_taxa_mrca_node
			all_taxa_mrca_node.remove_child(first_cluster_mrca_node)
			first_cluster_mrca_node.parent_node = None
			all_taxa_mrca_node.remove_child(second_cluster_leaf_node)
			second_cluster_leaf_node.parent_node = None
			# add these individual clusters' MRCA node as its children
			newnode.add_child(first_cluster_mrca_node)
			first_cluster_mrca_node.parent_node = newnode
			newnode.add_child(second_cluster_leaf_node)
			second_cluster_leaf_node.parent_node = newnode
			# update splits of the resulting tree
			Star_Tree_Initial.update_splits(delete_outdegree_one=False)
		
		else:
			"""
			case 4 - when both child clusters are leaf nodes 
			"""
			first_cluster_leaf_node = Star_Tree_Initial.find_node_with_taxon_label(clust_species_list[min_idx_i][0])
			second_cluster_leaf_node = Star_Tree_Initial.find_node_with_taxon_label(clust_species_list[min_idx_j][0])
			all_taxa_mrca_node = Star_Tree_Initial.mrca(taxon_labels=taxa_list)
			if (DEBUG_LEVEL >= 2):
				fp = open(Output_Text_File, 'a')
				fp.write('\n first cluster is a leaf - its label: ' + str(Node_Label(first_cluster_leaf_node)))      
				fp.write('\n second cluster is a leaf - its label: ' + str(Node_Label(second_cluster_leaf_node)))
				fp.write('\n label of all_taxa_mrca_node: ' + str(Node_Label(all_taxa_mrca_node)))
				fp.close()
			
			# create new internal node 
			newnode = dendropy.Node()
			# its parent node will be the previous MRCA node of all the taxa in two clusters
			all_taxa_mrca_node.add_child(newnode)
			newnode.parent_node = all_taxa_mrca_node
			all_taxa_mrca_node.remove_child(first_cluster_leaf_node)
			first_cluster_leaf_node.parent_node = None
			all_taxa_mrca_node.remove_child(second_cluster_leaf_node)
			second_cluster_leaf_node.parent_node = None
			# add these individual clusters' MRCA node as its children
			newnode.add_child(first_cluster_leaf_node)
			first_cluster_leaf_node.parent_node = newnode
			newnode.add_child(second_cluster_leaf_node)
			second_cluster_leaf_node.parent_node = newnode
			# update splits of the resulting tree
			Star_Tree_Initial.update_splits(delete_outdegree_one=False)

		#---------------------------------------------------------------------
		if (DEBUG_LEVEL >= 2):
			fp = open(Output_Text_File, 'a')
			fp.write('\n label of newnode: ' + str(Node_Label(newnode)))
			fp.write('\n label of all taxa mrca node (recomputed): ' + str(Node_Label(all_taxa_mrca_node)))      
			#fp.write('\n before inserting row col, Dist_Mat_clust_pair_NJ dimension: ' + str(Dist_Mat_clust_pair_NJ.size))
			fp.close()

		"""
		adjust the Dist_Mat_clust_pair_NJ by inserting one new row and column corresponding to the new cluster
		and then deleting the information of earlier two clusters
		"""
		# first append one row
		Dist_Mat_clust_pair_NJ = numpy.vstack((Dist_Mat_clust_pair_NJ, numpy.zeros((1, no_of_taxa_clust), dtype=numpy.float)))
		# then append one column
		Dist_Mat_clust_pair_NJ = numpy.hstack((Dist_Mat_clust_pair_NJ, numpy.zeros((no_of_taxa_clust + 1, 1), dtype=numpy.float)))
		# now reshape the distance matrix
		Dist_Mat_clust_pair_NJ = numpy.reshape(Dist_Mat_clust_pair_NJ, ((no_of_taxa_clust + 1), (no_of_taxa_clust + 1)), order='C')
		
		"""
		add taxa_list as a new element of clust_species_list
		"""
		clust_species_list.append(taxa_list)          
		
		"""
		now recompute the entries of this new row and column (which is indexed by no_of_taxa_clust), according to the NJ principle
		compute Dist_Mat_clust_pair_NJ[no_of_taxa_clust][m] entries where m != min_idx_i and m != min_idx_j
		"""
		for m in range(no_of_taxa_clust):
			if (m == min_idx_i) or (m == min_idx_j):
				continue
			Dist_Mat_clust_pair_NJ[no_of_taxa_clust][m] = (Dist_Mat_clust_pair_NJ[min_idx_i][m] + \
				Dist_Mat_clust_pair_NJ[min_idx_j][m] - Dist_Mat_clust_pair_NJ[min_idx_i][min_idx_j]) / 2.0
			Dist_Mat_clust_pair_NJ[m][no_of_taxa_clust] = Dist_Mat_clust_pair_NJ[no_of_taxa_clust][m]
		
		"""
		now remove the rows and columns corresponding to min_idx_i and min_idx_j
		"""
		Dist_Mat_clust_pair_NJ = numpy.delete(Dist_Mat_clust_pair_NJ, (min_idx_i), axis=0)	# delete the row
		Dist_Mat_clust_pair_NJ = numpy.delete(Dist_Mat_clust_pair_NJ, (min_idx_i), axis=1)	# delete the column
		Dist_Mat_clust_pair_NJ = numpy.delete(Dist_Mat_clust_pair_NJ, (min_idx_j - 1), axis=0)	# delete the row
		Dist_Mat_clust_pair_NJ = numpy.delete(Dist_Mat_clust_pair_NJ, (min_idx_j - 1), axis=1)	# delete the column

		# clear Norm_DistMat_ClustPair_NJ
		Norm_DistMat_ClustPair_NJ = numpy.delete(Norm_DistMat_ClustPair_NJ, (min_idx_i), axis=0)	# delete the row
		Norm_DistMat_ClustPair_NJ = numpy.delete(Norm_DistMat_ClustPair_NJ, (min_idx_i), axis=1)	# delete the column
		Norm_DistMat_ClustPair_NJ.fill(0)
		
		# remove individual clusters' taxa information from the clust_species_list
		clust_species_list.pop(min_idx_i)
		clust_species_list.pop(min_idx_j - 1)
		
		# decrement the number of clusters considered
		no_of_taxa_clust = no_of_taxa_clust - 1

	return
