# STAR_GLASS_STEAC
Implements species tree construction from incongruent gene trees with Incomplete Lineage Sorting (ILS) using a variety of existing methods - 1) STAR 2) STEAC (Liu et. al. 2009), 3) GLASS (Roch et. al. 2010)

Description
------------

Input
-----

A collection of gene trees with overlapping taxa set (sampled genes), having topological incongruence due to Incomplete Lineage Sorting (ILS).

Gene trees may be weighted (having branch length / coalescence time information) or may not be.

Output
-------

A species tree covering all the taxa of the gene trees. The species tree is aimed to be topologically closer to the model species tree (if available), or to the input gene trees.

Three methods are implemented and integrated in this framework, as follows:

1) STAR (Liu et. al. 2009) uses the coalescence rank of individual couplets to construct the species tree using a neighbor joining (NJ) approach.

2) STEAC (Liu et. al. 2009) uses the average coalescence time of individual couplets, computed with respect to the input gene trees. The coalescence time is approximated with the branch length information provided in the weighted input gene trees. This algorithm does not apply for unweighted gene trees.

3) GLASS (Roch et. al. 2010) uses the minimum coalescence time of individual couplets, computed with respect to the input gene trees. The coalescence time is approximated with the branch length information provided in the weighted input gene trees. This algorithm does not apply for unweighted gene trees.

Dependencies / Installation Requirements
-----------------------------------------

CSTBL is developed in Linux Systems (Ubuntu 14.04), using Python 2.7.

User needs to install following before using this package:

1) Python 2.7 (available in Ubuntu, by default)

Note: We have not tested the code on Python 3. Any user having Python 3 environment need to check the correct execution of our code, and optionally needs to upgrade it accordingly.

We plan to support Python 3 environment in some future release.

2) Dendropy 3.12.0 ( available on the link: https://pythonhosted.org/DendroPy/ )

Note: there is a new release of Dendropy 4.0 but we have used 3.12.0 for the implementation. We did not upgrade the code for Dendropy 4.0 support, so any user having this new version of Dendropy might need to check the functionalities of COSPEDBTree and possibly upgrade / replace / edit few dendrop[y related functions. So, we recommend users to use the earlier version of Dendropy, to avoid any conflict.

Support for Dendropy 4 and corresponding update of code will be done in a future release.

3) Numpy ( available on the link: http://www.numpy.org/ )

User can install Numpy using pip (python software downloader tool) module, which contains the latest Numpy module in it. We found that Numpy module in the traditional Apt-Get repository is of lower version.

UBUNTU version issues
----------------------

For systems having Ubuntu with lower versions (lower than 14.04), please notify in case of any errors due to OS environments.

Note: We do not support development version corresponding to Windows XP and MacOS, although that will be done in some future release.

Execution
------------

The file STAR_GLASS_STEAC.py is the main file in this package. It is to be executed with the following command line options, from a terminal. In terminal, go to the directory containing the source codes, and type the following commands:

chmod +x STAR_GLASS_STEAC.py (To change its permission to make it an executable file)

./STAR_GLASS_STEAC.py [options]

Details of the options are mentioned below:

-h, --help show this help message and exit

  -I INP_FILENAME, --INPFILE=INP_FILENAME
                        
                        name of the input file containing gene trees
  
  -O OUT_FILENAME, --OUTFILE=OUT_FILENAME
  
                        Name of the output file to contain target species tree. 
                        It is not mandatory to provide this option. If this option is not used, 
                        a folder named (STAR / GLASS / STEAC) will be created (depending upon the 
                        species tree generation method employed) in the directory containing the input 
                        gene tree list file. Within that new created directory, a file 'outtree_Newick.tre' will 
                        contain the output species tree.
                        
  -p INP_FILE_FORMAT, --inpform=INP_FILE_FORMAT
                        
                        1 - input file format is NEWICK (default)
                        2 - input file format is NEXUS
                        
  -m METHOD_USED, --method=METHOD_USED
                        
                        1 - average rank of the couplets (with respect to the
                        gene trees) will be used for NJ based clustering (STAR)                     
                        2 - minimum coalescence time of the couplets (with respect to the gene trees)
                        will be used for NJ based clustering (GLASS)
                        3 - average coalescence time of the couplets (with
                        respect to the gene trees) will be used for NJ based clustering (STEAC)

Example of a command (followed for the results published in the manuscript)

./STAR_GLASS_STEAC.py -I source_tree_input.txt -p1 -T supertree_topology_file.txt -m1 > out.txt

command descriptions:

1) -I specifies the input filename

2) source_tree_input.txt : contains the input collection o gene trees

3) -p option is for specifying the input tree format input file contains the trees in NEWICK format, as specified by the option (-p1) (1 stands for newick)

4) -m option is used to specify the species tree construction method. Here 1 indicates that STAR method is employed.

The output texts are printed at console. User can redirect the output results to any standard text file by using standard redirection operation (>). For example, in the above command, all the detailed results (textual descriptions) are redirected to file out.txt.

As mentioned in the command options, specification of output species tree containing file is not mandatory. In such a case, a folder named (STAR / GLASS / STEAC) will be created (depending upon the species tree generation method employed) in the directory containing the input gene tree list file. Within that new created directory, a file 'outtree_Newick.tre' will contain the output species tree.

Utilities
---------

All three species tree construction methods associate O(MN^2 + N^3) time complexity and O(N^2) space complexity, for N input taxa and M input trees.

For any queries, please contact
---------------

Sourya Bhattacharyya

Department of Computer Science and Engineering

Indian Institute of Technology Kharagpur

email: sourya.bhatta@gmail.com

