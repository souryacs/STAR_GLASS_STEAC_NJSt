#!/usr/bin/env python

##---------------------------------------------
''' 
this program is used to generate a spcies tree from a set of gene trees
the gene trees generally associate conflicts among the constituent genes (representing individual taxa)
in terms of the topological relationships
species tree follows simple approach to produce species tree with high performance index

Author: Sourya Bhattacharyya
Dept of CSE, IIT Kharagpur
V1.0 - 12.01.2015 - basic code
V2.0 - 19.06.2015 - removed Phylip NJ and metrics 4 and 5 implementations
''' 

## Copyright 2015 Sourya Bhattacharyya and Jayanta Mukherjee.
## All rights reserved.
##
## See "LICENSE.txt" for terms and conditions of usage.
##
##---------------------------------------------

import Header
from Header import *
import UtilFunc
from UtilFunc import *
import NJ_SpecTree
from NJ_SpecTree import *

##-----------------------------------------------------
# this function is useful to parse various options for input data processing
def parse_options():  
  parser = OptionParser()
    
  parser.add_option("-I", "--INPFILE", \
			  type="string", \
			  action="store", \
			  dest="INP_FILENAME", \
			  default="", \
			  help="name of the input file containing gene trees")
			  
  parser.add_option("-O", "--OUTFILE", \
			  type="string", \
			  action="store", \
			  dest="OUT_FILENAME", \
			  default="", \
			  help="name of the output file to contain target species tree")  
  
  parser.add_option("-p", "--inpform", \
			  type="int", \
			  action="store", \
			  dest="inp_file_format", \
			  default=1, \
			  help="1 - input file format is NEWICK (default) \
			  2 - input file format is NEXUS")
    
  parser.add_option("-m", "--method", \
			  type="int", \
			  action="store", \
			  dest="method_used", \
			  default=1, \
			  help="1 - average rank of the couplets (with respect to the gene trees) will be used for NJ based clustering (STAR) \
			  2 - minimum coalescence time of the couplets (with respect to the gene trees) will be used for NJ based clustering (GLASS) \
			  3 - average coalescence time of the couplets (with respect to the gene trees) will be used for NJ based clustering (STEAC)")
       
  opts, args = parser.parse_args()
  return opts, args
  
##-----------------------------------------------------
''' main function '''
def main():  
  opts, args = parse_options()

  ROOTED_TREE = True	#False 
  PRESERVE_UNDERSCORE = True 
  INPUT_FILE_FORMAT = opts.inp_file_format
  INPUT_FILENAME = opts.INP_FILENAME
  OUTPUT_FILENAME = opts.OUT_FILENAME
  METHOD_USED = opts.method_used
  
  global Output_Text_File
  
  if (INPUT_FILENAME == ""):
    print '******** THERE IS NO INPUT FILE SPECIFIED - RETURN **********'
    return
  else:
    print 'input filename: ', INPUT_FILENAME
    
  # according to the location of input filename
  # adjust the locations of the output files as well
  k = INPUT_FILENAME.rfind("/")
  if (k == -1):
    dir_of_inp_file = './'
  else:
    dir_of_inp_file = INPUT_FILENAME[:(k+1)]
  print 'dir_of_inp_file: ', dir_of_inp_file  
      
  if (OUTPUT_FILENAME == ""):
    # derive the output directory which will contain different output text results
    if (METHOD_USED == STAR):
      dir_of_curr_exec = dir_of_inp_file + 'STAR'
    elif (METHOD_USED == GLASS):
      dir_of_curr_exec = dir_of_inp_file + 'GLASS'
    elif (METHOD_USED == STEAC):
      dir_of_curr_exec = dir_of_inp_file + 'STEAC'      
    
    # append the current output directory in the text file
    Output_Text_File = dir_of_curr_exec + '/' + 'Complete_Desription.txt'
    # create the directory
    if (os.path.isdir(dir_of_curr_exec) == False):
      mkdr_cmd = 'mkdir ' + dir_of_curr_exec
      os.system(mkdr_cmd)               
    print 'Output_Text_File: ', Output_Text_File      
  else:
    # when we have specified the output file then corresponding directory becomes the dir_of_curr_exec
    k1 = OUTPUT_FILENAME.rfind("/")
    if (k1 == -1):
      dir_of_curr_exec = './'
    else:
      dir_of_curr_exec = OUTPUT_FILENAME[:(k1+1)]
    Output_Text_File = OUTPUT_FILENAME + '_complete_text_description'   
    # create the directory
    if (os.path.isdir(dir_of_curr_exec) == False):
      mkdr_cmd = 'mkdir ' + dir_of_curr_exec
      os.system(mkdr_cmd)                   
    print 'Output_Text_File: ', Output_Text_File
    
  fp = open(Output_Text_File, 'w')    
  
  fp.write('\n ================ status of options ================= (1 means ON)')
  fp.write('\n ROOTED_TREE: ' + str(ROOTED_TREE))
  fp.write('\n PRESERVE_UNDERSCORE: ' + str(PRESERVE_UNDERSCORE))
  fp.write('\n ===>>>  processing the file now ======== ')

  # note the program beginning time 
  start_timestamp = time.time()

  #-------------------------------------------------------------
  """ read the source trees collection and store it in a tree collection structure
  individual elements of this collection is thus a source tree """
  Gene_TreeList = Read_Gene_Data_Treelist(ROOTED_TREE, PRESERVE_UNDERSCORE, INPUT_FILE_FORMAT, INPUT_FILENAME)
  
  # these variables store the respective parameters of input gene trees
  tree_count = 0
  COMPLETE_INPUT_TAXA_LIST = []

  # first collect the taxa information from the source trees
  for tr_idx in range(len(Gene_TreeList)):
    tree_count = tree_count + 1
    # add - sourya
    # if this is an unweighted tree list 
    # then we cannot use coalescence time based processing 
    # and corresponding NJ based derivation    
    if (tree_count == 1):
      length_edges = Gene_TreeList[tr_idx].length()
      print '*** length edges: ', length_edges
      if (length_edges == 0):
	UNWEIGHTED_TREE = True
      else:
	UNWEIGHTED_TREE = False
      if (UNWEIGHTED_TREE) and ((METHOD_USED == GLASS) or (METHOD_USED == STEAC)):
	print '*** unweighted gene trees cannot use coalescence time information for species tree derivation -- return'
	return
    # end add - sourya
    taxa_labels_curr_tree = Gene_TreeList[tr_idx].infer_taxa().labels()
    fp.write('\n Tree no : ' + str(tree_count) + ' no of leaf nodes: ' + str(len(taxa_labels_curr_tree)))
    # append taxa information
    for t in taxa_labels_curr_tree:
      if t not in COMPLETE_INPUT_TAXA_LIST:
	COMPLETE_INPUT_TAXA_LIST.append(t)
	if (DEBUG_LEVEL > 2):
	  fp.write('\n new taxa : ' + str(t))
    
  # now process individual trees to find the couplet relations of those trees
  for tr_idx in range(len(Gene_TreeList)):
    DeriveCoupletRelations(Gene_TreeList[tr_idx], METHOD_USED)
    
  # now we print the information for all different couplets
  if (DEBUG_LEVEL >= 0):
    fp.write('\n  total no of taxa: ' + str(len(COMPLETE_INPUT_TAXA_LIST)))

  # closing the output text file before using it in a function call
  fp.close()

  ##-------------------------------------------
  ## check if the plotting of individual parameters are permitted
  #if (PLOTTING_ENABLED == 1):
    #for l in TaxaPair_Reln_Dict:
      #TaxaPair_Reln_Dict[l]._PlotStatistics(l, METHOD_USED, dir_of_curr_exec)

  #-------------------------------------------

  if 0:	#(DEBUG_LEVEL >= 2):
    fp = open(Output_Text_File, 'a')    
    fp.write('\n\n **** printing the couplet information ***')
    fp.close()
    for l in TaxaPair_Reln_Dict:
      TaxaPair_Reln_Dict[l]._PrintTaxaPairRelnInfo(l, METHOD_USED, Output_Text_File)  
      
  #-------------------------------------------------------------
      
  # generate a star network from the input taxa labels
  # form a newick formatted string containing the tree
  star_net_str = ""
  for i in range(len(COMPLETE_INPUT_TAXA_LIST)):
    star_net_str = star_net_str + str(COMPLETE_INPUT_TAXA_LIST[i])
    if (i < (len(COMPLETE_INPUT_TAXA_LIST) - 1)):
      star_net_str = star_net_str + ","
  star_net_str = "(" + star_net_str + ")"

  # this tree denotes the initial star configuration (rooted at the central hub)
  Star_Tree_Initial = dendropy.Tree.get_from_string(star_net_str, schema="newick", \
						    preserve_underscores=PRESERVE_UNDERSCORE, \
						    default_as_rooted=True)          
  
  fp = open(Output_Text_File, 'a')    
  fp.write('star_net_str: ' + str(star_net_str))
  fp.write('\n from tree ---: ' + Star_Tree_Initial.as_newick_string())
  fp.close()
  
  # now perform the agglomerative clustering technique based on the extra lineages
  Form_Species_Tree_NJ_Cluster(Star_Tree_Initial, COMPLETE_INPUT_TAXA_LIST, METHOD_USED, Output_Text_File)
  
  # note the time
  end_timestamp = time.time()      
  
  fp = open(Output_Text_File, 'a')
  fp.write('\n --- output species tree (in newick format): ' + Star_Tree_Initial.as_newick_string())    
  fp.close()
  
  # write this tree on a separate text file
  if (OUTPUT_FILENAME == ""):
    out_treefilename = dir_of_curr_exec + '/' + 'outtree_Newick.tre'
  else:
    out_treefilename = OUTPUT_FILENAME
  
  # we write the unweighted supertree
  outfile = open(out_treefilename, 'w')
  outfile.write(Star_Tree_Initial.as_newick_string())	
  outfile.close()
  
  # we write the time associated with the execution of this method
  fp = open(Output_Text_File, 'a')
  fp.write('\n \n\n ===============>>>>>>>>>>>>>>> TIME COMPLEXITY OF THE METHOD (in seconds) ')
  fp.write('\n \n Total time taken (in seconds) : ' + str(end_timestamp - start_timestamp))
  fp.close()
  
#-----------------------------------------------------
if __name__ == "__main__":
    main() 




