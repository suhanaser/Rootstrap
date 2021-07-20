# rootstrap.py: A script to calculate the rootstrap support values for all the branches
A script to calculate rootstrap values in Python. Please note that it's probably more convenient to do this in IQ-TREE: http://www.iqtree.org/doc/Rootstrap
This script was written before the IQ-TREE implementation, and is what we used in the paper. It also serves as a good check on whether the IQ-TREE implementation works. Both implementations should give you the same answer (they do for us in all of our tests!).
### Installation
To run this script you need Python 3 and three dependencies: `ete3`, `Biopython`, and `dendropy`. You can install these however you like, but if you use Anaconda, you can create an environment and install them as follows:
```
conda create -n rootstrap python=3.6
conda activate rootstrap
conda install -c etetoolkit ete3
conda install -c conda-forge biopython
conda install -c bioconda dendropy
```
### Running the script
You can either download just the `rootstrap.py` script from this repo, or you can get the whole repo. It's proabably easiest to do the latter, like this:
```
mkdir rootstrap
cd rootstrap
git clone git clone https://github.com/suhanaser/Rootstrap.git
```
Now you can run the script as follows
```
python Rootstrap/rootstrap.py
```
To get the help information, just do this:
```
python Rootstrap/rootstrap.py -h
```
That should show you the following help:
```
Syntax:
rootstrap.py <tree file> <bootstrap trees file> <is rooted> <outgroup file>
<treefile>  		 The tree file where you want to calculate the rootstrap support values in Newick format (e.g. tree.treefile)
.
<bootstrap trees file>	 All the bootstrap trees in Newick format. e.g. .ufboot file from IQ-TREE (tree.ufboot)
<is rooted>  		 Are the ML tree and the boostrap trees rooted or not. default: True
			 True - The trees assumed to be rooted and no outgroup taxa infromation is required
			 False - The trees assumed to be unrooted and outgroup taxa infromation is required
<outgroup file> 	 outgroup taxa in Nexus format.
			 The outgroup block can be part of the alignment file or in a separate file
```
An example of how to use the script is as follows. This uses test files contained in the repository itself, so to run this you'll need to have downloaded the full repository as above.
```
python Rootstrap/rootstrap.py Rootstrap/tree.treefile Rootstrap/tree.ufboot
```

This analysis should produce a file called `tree.rootstrap`. This is a tree in nexus format, where each node is now labelled with the rootstrap proportion e.g. `[&rootstrap=0.1]`. You can open this file in FigTree and look at the rootstrap percentages. 

### rBED_rSED.py: A script to calculate the root Branchlength Error Distance (rBED) and the root Split Error Distance (rSED)
Syntax:
`rBED_rSED.py <tree file>`

`<treefile>`: The tree file (rooted tree only in Newick format) where you want to calculate rBED and rSED values.

Note: In order to calculate rBED and rSED values, the true root should be known (or assumed to be known) in advance.
