# STAR_GLASS_STEAC
Implements species tree construction from incongruent gene trees with Incomplete Lineage Sorting (ILS). 

Following three different methods are implemented in this package.

1) STAR (Liu et. al. 2009), which constructs species trees using couplet coalescence rank information.

2) STEAC (Liu et. al. 2009), which constructs species trees using average coalescence time for individual couplets.

3) GLASS (Roch et. al. 2010), which constructs species trees using minimum coalescence time for individual couplets.

4) NJ_st (Liu et. al. 2011), which constructs species trees using average branch count for individual couplets.

All of these methods use NJ based tree construction technique. Distance matrix (at couplet level) is 
constructed using above mentioned measures.

-----------------
Description
----------------

Input
----------

A collection of gene trees with overlapping taxa set (sampled genes), having topological incongruence due to Incomplete Lineage Sorting (ILS).

Gene trees may be weighted (having branch length / coalescence time information) or may not be.

Support for multi allele gene trees
---------------------------------------------

From version 2.0, the code supports multi allele gene trees. 

Note: Basic dendropy supports single allele gene trees. 
So, we have introduced a parser algorithm which takes care of the alleles and converts input multi 
allele gene trees into a single allele one, suitable for processing in dendropy.

To support multi allele gene trees, we support only newick formatted input gene trees for the moment.
So, user needs to provide newick format gene tree list files as input. Any other format is not 
supported for the moment.

Output
-----------

A species tree covering all the taxa of the gene trees. The species tree is aimed to be topologically closer to the 
model species tree (if available), or to the input gene trees.

Methods implemented
-------------------

Three methods are implemented and integrated in this framework, as follows:

1) STAR (Liu et. al. 2009) uses the coalescence rank of individual couplets to construct the 
species tree using a neighbor joining (NJ) approach.

2) STEAC (Liu et. al. 2009) uses the average coalescence time of individual couplets, computed with respect to the input gene trees. 
The coalescence time is approximated with the branch length information provided in the weighted input gene trees. 
This algorithm does not apply for unweighted gene trees.

3) GLASS (Roch et. al. 2010) uses the minimum coalescence time of individual couplets, computed 
with respect to the input gene trees. The coalescence time is approximated with the branch length information provided 
in the weighted input gene trees. This algorithm does not apply for unweighted gene trees.

Features
----------
1) Couplet based processing. All three species tree construction methods associate O(MN^2 + N^3) time complexity 
and O(N^2) space complexity, for N input taxa and M input trees.

2) Single allelle gene trees are provided as inputs.


Dependencies / Installation Requirements
-----------------------------------------

This package is developed in Linux Systems (Ubuntu 14.04), using Python 2.7.

User needs to install following before using this package:

1) Python 2.7 (available in Ubuntu, by default)

Note: We have not tested the code on Python 3. Any user having Python 3 environment need to check the correct 
execution of our code, and optionally needs to upgrade it accordingly.

We plan to support Python 3 environment in some future release.

2) Dendropy 3.12.0 ( available on the link: https://pythonhosted.org/DendroPy/ )

**** Note: there is a new release of Dendropy 4.0 but we have used 3.12.0 for the implementation. We did not upgrade 
the code for Dendropy 4.0 support, so any user having this new version of Dendropy might need to check the 
functionalities of this package and possibly upgrade / replace / edit few dendrop[y related functions. So, we 
recommend users to use the earlier version of Dendropy, to avoid any conflict.

Support for Dendropy 4 and corresponding update of code will be done in a future release.

3) Numpy ( available on the link: http://www.numpy.org/ )

User can install Numpy using pip (python software downloader tool) module, which contains the latest Numpy module in it. 
We found that Numpy module in the traditional Apt-Get repository is of lower version.

UBUNTU version issues
----------------------

For systems having Ubuntu with lower versions (lower than 14.04), please notify in case of any errors due to OS environments.

Note: We do not support development version corresponding to Windows XP and MacOS, although that will be done in some future release.

Execution
------------

The file STAR_GLASS_STEAC.py is the main file in this package. It is to be executed with the following 
command line options, from a terminal. In terminal, go to the directory containing the source codes, and 
type the following commands:

chmod +x STAR_GLASS_STEAC.py (To change its permission to make it an executable file)

./STAR_GLASS_STEAC.py [options]

Details of the options are mentioned below:

-h, --help show this help message and exit

  -I INP_FILENAME, --INPFILE=INP_FILENAME
                        
                        name of the input file containing gene trees
  
  -O OUT_FILENAME, --OUTFILE=OUT_FILENAME
  
                        Name of the output file to contain target species tree. 
                        It is not mandatory to provide this option. If this option is not used, 
                        a folder named (STAR / GLASS / STEAC/NJ_st) will be created (depending upon the 
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
                        4 - average branch count of the couplets (with
                        respect to the gene trees) will be used for NJ based clustering (NJ_st)

Example of a command (followed for the results published in the manuscript)

./STAR_GLASS_STEAC.py -I source_tree_input.txt -p1 -T supertree_topology_file.txt -m1 -g2 > out.txt

command descriptions:

1) -I specifies the input filename

2) source_tree_input.txt : contains the input collection of gene trees

3) -p option is for specifying the input tree format input file contains the trees in NEWICK format, 
as specified by the option (-p1) (1 stands for newick)

Note: For multi allele gene trees, p should be provided as 1.

4) -m option is used to specify the species tree construction method. Here 1 indicates that STAR method is employed.

5) -g option is used to specify whether the gene trees contain multiple allele information. 
In such a case, -g2  option is used.

The output texts are printed at console. User can redirect the output results to any standard text file by using 
standard redirection operation (>). For example, in the above command, all the detailed results (textual descriptions) 
are redirected to file out.txt.

As mentioned in the command options, specification of output species tree containing file is not mandatory. 
In such a case, a folder named (STAR / GLASS / STEAC/NJ_st) will be created (depending upon the species tree 
generation method employed) in the directory containing the input gene tree list file. Within that new created 
directory, a file 'outtree_Newick.tre' will contain the output species tree. Another file named 
'Complete_Description.txt' will be created within that directory. It will contain the details of execution.

For any queries, please contact
----------------------------------------

Sourya Bhattacharyya

Department of Computer Science and Engineering

Indian Institute of Technology Kharagpur

email: sourya.bhatta@gmail.com


