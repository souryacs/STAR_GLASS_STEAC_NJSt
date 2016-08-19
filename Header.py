#!/usr/bin/env python

"""
program for species tree estimation using features computed for individual couplets
computes four criteria - 
1) Coalescence Rank statistics (STAR)
2) Min coalescence time (GLASS)
3) Avg coalescence time (STEAC)
4) Internode count statistics (NJst)
"""

import dendropy
from optparse import OptionParser
import math
import time
import os
import numpy
#import matplotlib.pyplot as plt
#import pylab
#import re

""" 
this dictionary contains the features for individual couplets
"""
TaxaPair_Reln_Dict = dict()

""" 
this list contains the complete set of taxa present in the input trees 
"""
COMPLETE_INPUT_TAXA_LIST = []

"""
this is the allelle dictionary corresponding to the taxon names and the allelles
"""
Allele_Dict = dict()

"""
this is a dictionary indexed by tree number
individual trees contain those taxon information which contain multiple allele 
a tree contain a list, where individual element of a list is a sublist of 2 elements:
1) taxon of multiple allele information
2) count of multiple allele
"""
Tree_Allele_Dict = dict()

"""
this text file stores all the printing output
"""
Output_Text_File = 'complete_output_description.txt'

"""
this is the debug level
set for printing the necessary information
"""
DEBUG_LEVEL = 0

## this is for plotting various statistical parameters
## computed using the input gene trees
#PLOTTING_ENABLED = 0

# methods implemented in the current package
STAR = 1
GLASS = 2
STEAC = 3
NJ_st = 4

# type of input gene trees
SINGLE_ALLELE = 1
MULTI_ALLELE = 2

#-----------------------------------------------------
""" 
this class defines the connectivity relationship between a pair of taxa
initially the information are obtained from the input source trees
"""
class Reln_TaxaPair(object):  
	def __init__(self):
		"""
		for single allelle gene trees, 
		this is the count of trees for which the couplet is supported
		for multi allelle gene trees, this is the number of times this couplet occurs for all the gene trees
		here, one gene tree may contain multiple instances of the same couplet
		"""
		self.tree_support_count = 0        
		""" 
		for single allelle gene trees, this list corresponds to the LCA rank of this couplet for individual source trees
		for multi allelle gene trees, a particular gene tree may contain more than one instance of this couplet
		individual LCA rank values, in such a case, are inserted in this list
		"""
		self.rank_statistics_list = []
		"""
		branch length distance between individual leaves to their MRCA nodes
		for single allelle gene trees, this list corresponds to the coalescence time information of this couplet for individual source trees
		for multi allelle gene trees, a particular gene tree may contain more than one instance of this couplet
		individual coalescence information, in such a case, are inserted in this list    
		"""
		self.coalescence_time_list = []
		""" 
		This variable contains the sum of internode count distances between this couplet
		with respect to all the input gene trees
		"""
		self.level_info_list = []
		
	## histogram plotting
	#def PlotHist(self, inp_key, inp_no, outdir, inp_str):
		#title_str = 'Histogram of ' + inp_str + ' for the couplet ' + str(inp_key[0]) + ' and ' + str(inp_key[1])
		#if (inp_no == 11):
			#arr = self.rank_statistics_list
		#elif (inp_no == 21) or (inp_no == 31):
			#arr = self.coalescence_time_list
		#elif (inp_no == 41):
			#arr = self.extra_lineage_sum_list
		#elif (inp_no == 42):
			#arr = self.lineage_sum_list
		#elif (inp_no == 43):
			#arr = self.level_info_list
		#elif (inp_no == 44):
			#arr = self.score_product_xl_level_list
		#elif (inp_no == 51):
			#arr = self.accumulated_rank_list
			
		## the histogram of the data
		#n, bins, patches = plt.hist(arr, bins=100, normed=0, facecolor='g', alpha=0.75, hold=False)
		#plt.xlabel(inp_str)
		#plt.ylabel('Count')
		#plt.title(title_str)
		##plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
		##plt.axis([40, 160, 0, 0.03])
		#plt.grid(True)
		##plt.show()  
		#output_image_filename = outdir + inp_str + '_' + str(inp_key[0]) + '_and_' + str(inp_key[1]) + '.jpg'
		#pylab.savefig(output_image_filename)	# saves the figure in custom specified image
	
	"""
	this function adds the count of tree if that tree supports the current couplet 
	under consideration
	"""
	def _IncrSupportTreeCount(self):
		self.tree_support_count = self.tree_support_count + 1
	
	"""
	this function returns the number of trees supporting the couplet
	"""
	def _GetSupportTreeCount(self):
		return self.tree_support_count        
	
	"""
	adds the internode count distance for this couplet, with respect to a particular input tree
	"""
	def _AddLevel(self, val):
		self.level_info_list.append(val)

	"""
	returns the average internode count distance between this couplet
	"""
	def _GetAvgSumLevel(self):
		return (sum(self.level_info_list) * 1.0) / self.tree_support_count
	
	"""
	this function adds the LCA based coalescence rank information of this taxa pair
	corresponding to one input gene tree
	"""
	def _AddRankInfo(self, rank_val):
		self.rank_statistics_list.append(rank_val)
	
	"""
	computes the average coalescence rank information for this taxa pair
	used for the STAR method
	"""
	def _GetAvgRank(self):
		return (sum(self.rank_statistics_list) * 1.0) / self.tree_support_count

	"""
	Appends the coalescence time information for this couplet
	with respect to a particular gene tree
	the time value (coalescence_sum) is computed for both the taxa, with respect to their LCA node in the current gene tree
	"""
	def _AddCoalescenceInfo(self, coalescence_sum):
		self.coalescence_time_list.append(coalescence_sum)  

	"""
	computes the average coalescence time for this taxa pair
	used for the STEAC method
	"""
	def _GetAvgCoalescenceTime(self):
		return (sum(self.coalescence_time_list) * 1.0) / len(self.coalescence_time_list)
	
	"""
	computes the minimum coalescence time information for this taxa pair  
	used for the GLASS method
	"""
	def _GetMinCoalescenceTime(self):
		return min(self.coalescence_time_list)

	"""
	this function prints information for the current couplet
	"""
	def _PrintTaxaPairRelnInfo(self, key, METHOD_USED, out_text_file):
		fp = open(out_text_file, 'a')    
		fp.write('\n Current couplet: ' + str(key))
		if (METHOD_USED == STAR):
			fp.write('\n *** avergae coalescence rank: ' + str(self._GetAvgRank()))
		elif (METHOD_USED == GLASS) or (METHOD_USED == STEAC):
			fp.write('\n *** coalescence time info --- min: ' + str(self._GetMinCoalescenceTime()) + ' avg: ' + str(self._GetAvgCoalescenceTime()))
		elif (METHOD_USED == NJ_st):
			fp.write('\n *** Average internode count: ' + str(self._GetAvgSumLevel()))
		fp.close()
		
	## this function plots various statistics distribution
	#def _PlotStatistics(self, key, METHOD_USED, out_dir):
		#if (METHOD_USED == STAR):
			#self.PlotHist(key, 11, out_dir, 'coalescence rank')
		#if (METHOD_USED == GLASS): 
			#self.PlotHist(key, 21, out_dir, 'coalescence time')
		#if (METHOD_USED == STEAC):
			#self.PlotHist(key, 31, out_dir, 'coalescence time')
		#if (METHOD_USED == 4):
			#self.PlotHist(key, 41, out_dir, 'extra lineage sum')
			#self.PlotHist(key, 42, out_dir, 'total lineage sum')
			#self.PlotHist(key, 43, out_dir, 'level (branch) count')
			#self.PlotHist(key, 44, out_dir, 'product XL and level')
		#if (METHOD_USED == 5):
			#self.PlotHist(key, 51, out_dir, 'accumulated coalescence rank')
			
		#return
      