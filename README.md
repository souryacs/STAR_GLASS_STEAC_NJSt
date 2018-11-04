STAR_GLASS_STEAC
---------------

Implements various species tree construction methods from incongruent gene trees with Incomplete Lineage Sorting (ILS). Following four methods are implemented in this package.

1) STAR (Liu et. al. 2009), which constructs species trees using couplet coalescence rank information.

2) STEAC (Liu et. al. 2009), which constructs species trees using average coalescence time for individual couplets.

3) GLASS (Roch et. al. 2010), which constructs species trees using minimum coalescence time for individual couplets.

4) NJ_st (Liu et. al. 2011), which constructs species trees using average internode count for individual couplets.

All of these methods use NJ based tree construction technique. Distance matrix (at couplet level) is constructed using any of the above mentioned measures.


Description
----------------

Input
----------

A collection of gene trees with overlapping taxa set (sampled genes), having topological incongruence due to Incomplete Lineage Sorting (ILS). For execution of the methods STEAC or GLASS, the gene trees may be weighted (having branch length / coalescence time information). Other methods do not use weight information of the gene trees.

Support for multi allele gene trees
---------------------------------------------

From version 2.0, the code supports multi allele gene trees. 

Note: The code uses python package Dendropy (Sukumaran et. al. 2010) (version 3.12.0) for basic processing 
of phylogenetic trees. As the package supports only single allele gene trees, we have adapted the current package 
for processing multi-allele gene trees. Specifically, we have introduced a novel parser algorithm to detect the alleles, 
and to convert the input multi-allele gene trees into Dendropy compatible single allele trees.

*** Currently multi allele gene trees, if provided, should be given in newick format. Any other format (including nexus) 
is not supported. For single allele gene trees, both newick and nexus format are supported.

Output
-----------

A species tree covering all the taxa of the gene trees. 

Methods implemented
-------------------

Four methods are implemented and integrated in this framework, as follows:

1) STAR (Liu et. al. 2009) uses the coalescence rank of individual couplets to construct the species tree using a neighbor joining (NJ) approach.

2) STEAC (Liu et. al. 2009) uses the average coalescence time of individual couplets, computed with respect 
to the input gene trees. The coalescence time is approximated with the branch length information provided 
in the weighted input gene trees. This algorithm requires weighted gene trees.

3) GLASS (Roch et. al. 2010) uses the minimum coalescence time of individual couplets, computed with respect to the 
input gene trees. The coalescence time is approximated with the branch length information provided 
in the weighted input gene trees. This algorithm requires weighted gene trees.

4) NJ_st (Liu et. al. 2011): Here, the couplet based average internode distance measure is employed for NJ based 
species tree construction.

Complexity
----------

All four species tree construction methods associate O(MN^2 + N^3) time complexity and O(N^2) space complexity, 
for N input taxa and M input trees.

Dependencies / Installation Requirements
-----------------------------------------

This package is developed in Linux Systems (Ubuntu 14.04), using Python 2.7. User may download the zipped archieve from 
GitHub or clone the existing repository. User needs to also install following before using this package:

1) Python 2.7 (available in Ubuntu, by default)

Note: The current version does not support Python 3, and a future release is intended to support it.

2) Dendropy 3.12.0 ( available on the link: https://pythonhosted.org/DendroPy/ )

**** Note: there is a new release of Dendropy 4.1.0 but we have used 3.12.0 for the implementation. We did not upgrade 
the code for Dendropy 4.1.0 support, so user needs to use the older version of dendropy to use this code. 
Support for Dendropy 4 and corresponding update of code will be done in a future release.

3) Numpy ( available on the link: http://www.numpy.org/ )

User can install Numpy using pip (python software downloader tool) module, which contains the latest module.

UBUNTU version issues
----------------------

For systems having Ubuntu with lower versions (lower than 14.04), please notify in case of any errors due to 
OS environments.

Note: We do not support development version corresponding to Windows XP and MacOS, although that will be 
done in some future release.

Execution
------------

The file STAR_GLASS_STEAC.py is the main file in this package. It is to be executed with the following 
command line options:

python STAR_GLASS_STEAC.py [options]

Details of the options are mentioned below:

-h, --help 
	
	show this help message and exit

  -I INP_FILENAME, --INPFILE=INP_FILENAME
                        
	name of the input file containing gene trees
  
  -O OUT_FILENAME, --OUTFILE=OUT_FILENAME
  
	Name of the output file to contain target species tree. 
	It is not mandatory to provide this option. If this option is not used, 
	a folder named (STAR / GLASS / STEAC / NJ_st) will be created (depending upon the 
	species tree generation method employed) in the directory containing the input 
	gene tree list file. Within that new created directory, a file 'outtree_Newick.tre' will 
	contain the output species tree.
                    
  -g GENE_TYPE --genetype=GENE_TYPE
			  
	1 - gene trees are single allele type (default) 
	2 - gene trees are multi allele type
	User need to specify this option for multi allele gene trees.
                          
  -m METHOD_USED, --method=METHOD_USED
                        
	1 - average rank of the couplets (with respect to the
	gene trees) will be used for NJ based clustering (STAR)                     
	2 - minimum coalescence time of the couplets (with respect to the gene trees)
	will be used for NJ based clustering (GLASS)
	3 - average coalescence time of the couplets (with
	respect to the gene trees) will be used for NJ based clustering (STEAC)
	4 - average internode count of the couplets (with
	respect to the gene trees) will be used for NJ based clustering (NJ_st)

Example of a command (followed for the results published in the manuscript)

Example 1: Single allele gene trees with GLASS method

python STAR_GLASS_STEAC.py -I source_tree_input.txt -O out_species_tree.tre -g1 -m2 

Example 2: Multi allele gene trees with method NJst

python STAR_GLASS_STEAC.py -I source_tree_input.txt -O out_species_tree.tre -g2 -m4

*** Note: For multi allele gene trees, option -p should be set as 1 (newick format).

*** Note: if the output species tree is not provided as an input option, a folder  named (STAR / GLASS / STEAC / NJ_st) 
will be created (depending upon the species tree generation method) in the same directory containing the input gene 
tree list. A file 'outtree_Newick.tre' within the created directory would contain the output species tree. 
Another file named 'Complete_Description.txt' within this folder would contain the detailed log of execution.

Citation
--------

Upon using the packge, user should cite the following articles:

1) Sourya Bhattacharyya, Jayanta Mukherjee, IDXL: Species Tree Inference Using Internode Distance and Excess Gene Leaf Count, Journal of Molecular Evolution (Springer), volume 85, issue 1-2, pp. 57-78, 2017, DOI: 10.1007/s00239-017-9807-7

2) Sourya Bhattacharyya, Jayanta Mukhopadhyay, Accumulated Coalescence Rank and Excess Gene Count for Species Tree Inference, proceedings of 3rd International Conference on Algorithms for Computational Biology (AlCoB), Trujillo, Spain, Springer LNBI 9702, pp. 93-105.


For any queries, please contact
----------------------------------------

Sourya Bhattacharyya

La Jolla Institute for Allergy and Immunology

La Jolla, CA 92037, USA

email: sourya.bhatta@gmail.com


