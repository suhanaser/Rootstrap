rootstrap.py: A script to calculate the rootstrap support values for all the branches

Syntax:
rootstrap.py <tree file> <bootstrap trees file> <is rooted> <outgroup file>
<treefile>              : The tree file where you want to calculate the rootstrap support values in Newick format (e.g. tree.treefile).
                          The tree can be rooted or unrooted.
<bootstrap trees file>  : All the bootstrap trees in Newick format. 
                          For example: .ufboot file from IQ-TREE (tree.ufboot)
<is rooted>             : Are the ML tree and the boostrap trees rooted or not.
                          True - The trees assumed to be rooted and no outgroup taxa infromation is required
                          False - The trees assumed to be unrooted and outgroup taxa infromation is required
                          (default: True)
<outgroup file>         : outgroup taxa in Nexus format. 
                          The outgroup block can be part of the alignment (e.g. OG_File1.nex) file or in a separate file (e.g. OG_File2.nex)
