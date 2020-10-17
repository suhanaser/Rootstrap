`rootstrap.py`: A script to calculate the rootstrap support values for all the branches

Syntax:

```bash
rootstrap.py <tree file> <bootstrap trees file> <is rooted> <outgroup file>

<treefile>                     The tree file where you want to calculate the rootstrap support values in Newick format (e.g. tree.treefile). 
                               The tree can be rooted or unrooted.

<bootstrap trees file>         All the bootstrap trees in Newick format. 
                               For example: .ufboot file from IQ-TREE (tree.ufboot)

<is rooted>                    Are the ML tree and the boostrap trees rooted or not. 
                               True - The trees assumed to be rooted and no outgroup taxa infromation is required  
                               False - The trees assumed to be unrooted and outgroup taxa infromation is required  
                               (default: True)

<outgroup file>                outgroup taxa in Nexus format. 
                               The outgroup block can be part of the alignment (e.g. OG_File1.nex) file or in a separate file (e.g. OG_File2.nex)
```

-------------------------------------------------------------------------------------
`rBED_rSED.py`: A script to calculate the root Branchlength Error Distance (rBED) and the root Split Error Distance (rSED)

Syntax:

```bash
rBED_rSED.py                   <tree file>

<treefile>                     The tree file (rooted tree only in Newick format) where you want to calculate rBED and rSED values.
```

Note: In order to calculate rBED and rSED values, the true root should be known (or assumed to be known) in advance.

## Installation

For running the scripts in this repository you will need Python 3, and the dependencies
listed in `requirements.txt`.

It is recommended to use a virtual environment to avoid complications with system
dependencies and Python.

For users that prefer `pip`, the following should work assuming you have Python 3
`python` and `pip` in your command line `$PATH` (if using an older OS, you may still
have `python3` and `pip3`, so just use those commands instead).

```bash
$ python -m venv venv             # create a virtual environment called venv with pip
$ source venv/bin/activate        # load the virtual environment, note the change in some shells displaying the venv name
(venv) $ pip install -r requirements.txt
# ... it should take a few minutes to download dependencies from PYPI.org and install in your venv (not in your OS Python modules)
(venv) $ python rootstrap.py 
Error: Please provide the ML tree file and the bootstrap trees file
(venv) $ python rBED_rSED.py 
Error: Please provide the rooted ML tree file
```

If you had a similar output to the commands above, your installation worked successfully.
Have a read at the command syntax instructions to use the scripts.

For Conda users, follow their instructions for [managing environments](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)
and just run `pip`. Alternativelly, you can consult [Conda Forge](https://conda-forge.org/)
searching for equivalent packages and installing them via `conda install $package-name`.
