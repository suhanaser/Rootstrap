# -*- coding: utf-8 -*-
"""
Created on Sun Jul 19 15:35:48 2020

@author: Suha Naser-Khdour
"""

import sys
from ete3 import Tree
from Bio.Nexus import Nexus
from collections import OrderedDict
import itertools as ite
import os
import dendropy

def all_possible_roots(treefile):
    '''
    Parameters
    ----------
    treefile: The ML tree (in newick format)
    
    Returns
    -------
    roots: rooted ML tree with all possible root placements
    '''
    t = Tree(treefile)
    roots = OrderedDict()
    k = 1
    for n in t.traverse():
        if not n.is_root():
            if not n.is_leaf():
                n.add_features(name='n'+str(k))
                k += 1
        else:
            n.add_features(name='n0')
        
    for n in t.traverse():
        if not n.is_root():
            t.set_outgroup(n)
            roots[n.name] = t.write(format=9)
        else:
            roots[n.name] = t.write(format=9)
    return roots

def Read_Nex(f):
    '''
    Parameters
    ----------
    f : outgroup file in Nexus format

    Returns
    -------
    og : list of outgroup taxa

    '''
    og = []
    with open(f, 'r') as inf:
        for line in inf:
            if ('TAXSET outgroups' in line) or ('taxset outgroups' in line):
                og = line[line.index('=')+2:line.index(';')].split(" ")
    return og

def caluclate_rootstrap(treeFile, bootFile, is_rooted, out_group):
    '''
    Parameters
    ----------
    treeFile: rooted tree in newick format (.treefile in IQ-TREE)
    bootFile: rooted bootstrap trees in newick format (e.g. .ufboot file in IQ-TREE)
    rooted: if the bootstrap trees are rooted (defult is True). If not rooted provide outgroup taxa file
    og: A file with outgroup taxa in Nexus format
    
    Returns
    -------
    rootstrapTree: rooted tree with rootstrap support values as branch lengths in newick format
    '''

    boottrees = []
    trees = []
    polyphyly = 0
    N_boottrees = 0
    if not is_rooted:
        if out_group == None:
            raise SystemExit('Error: Please provide outgroup taxa in Nexus format')
        ML_tree = Tree(treeFile)
        try:
            og = Read_Nex(out_group) #get the outgroup taxa
        except:
            raise SystemExit('Error: Cannot find outgroup taxa')
        if len(og) == 1: #if there is one outgroup taxon use it to root the tree
            ML_root = ML_tree.search_nodes(name=og[0])[0]
        else: #if there are more than one outgroup taxon find their common ancestor
            ML_root = ML_tree.get_common_ancestor(og)
        if not ML_root.is_root():
            ML_tree.set_outgroup(ML_root)
        ingroup = [n.name for n in ML_tree.get_leaves() if n.name not in og]
        try:#check if the ingroup is monophyletic
            if ML_tree.check_monophyly(values=ingroup, target_attr="name", ignore_missing=True)[0]:
                ML_tree.prune(ingroup) #prune ingroup taxa only
                rootedMLtree = os.path.splitext(treeFile)[0]+'_rooted.treefile'
                ML_tree.write(outfile=rootedMLtree) #write the rooted ML tree with ingroup taxa only to a file
            else:
                 raise SystemExit('Error: ML ingroup taxa are not monophyletic')
        except:
                    raise SystemExit('Error: ML ingroup taxa are not monophyletic')
                
        with open(bootFile, 'r') as f:
            for tree in f:
                N_boottrees += 1
                t = Tree(tree)
                ingroup = [n.name for n in t.get_leaves() if n.name not in og]
                if len(og) == 1: #if there is one outgroup taxon use it to root the tree
                    root = t.search_nodes(name=og[0])[0]
                elif len(og) > 1: #if there are more than one outgroup taxon find their common ancestor
                    root = t.get_common_ancestor(og)
                else: #if there is no outgroup taxa raise an error
                    raise SystemExit('Error: Please provide outgroup taxa in Nexus format')
                if not root.is_root():
                    t.set_outgroup(root)
                try:#check if the ingroup is monophyletic
                    if t.check_monophyly(values=ingroup, target_attr="name", ignore_missing=True)[0]:
                        trees.append(t.write(format=9))
                    else:
                        polyphyly += 1
                except:
                    polyphyly += 1
        
        for tree in trees:
            t = Tree(tree)
            t.prune(ingroup)
            boottrees.append(t.write(format=9))
    else: #If you are using rooted ML tree and rooted bootstrap trees (e.g. NR model)
        ML_tree = Tree(treeFile)
        with open(bootFile, 'r') as f:
            for tree in f:
                N_boottrees += 1
                t = Tree(tree)
                boottrees.append(t.write(format=9))

    booted = [(g[0], len(list(g[1]))) for g in ite.groupby(boottrees)] #a list of all unique bootstrap trees with thier number of occurrence
    boottrees = []
    for b in booted:
        t2 = Tree(b[0])
        x = []
        for n in t2.traverse():
            if n.is_root():
                for child in n.children:
                    if child.is_leaf():
                        x.append([child.name])
                    else:
                        x.append([i.name for i in child.get_descendants()])
                boottrees.append([b[1],x])
    if is_rooted:
        roots = all_possible_roots(treeFile)
    else:
        roots = all_possible_roots(rootedMLtree)

    rootstrap_value = dict.fromkeys(roots.keys(), 0)
    for node, rooted in roots.items():
        t1 = Tree(rooted)
        x = []
        for n in t1.traverse():
            if n.is_root():
                for child in n.children:
                    if child.is_leaf():
                        x.append([child.name])
                    else:
                        x.append([i.name for i in child.get_descendants()])
        y = [set(i) for i in x]
        for split in boottrees:
            z = [set(i) for i in split[1]]
            if len(y) == len(z):
                for group in y:
                    if group in z:
                        z.remove(group)
            if len(z) == 0:
                rootstrap_value[node] += split[0]/N_boottrees
            else:
                rootstrap_value[node] += 0

    if is_rooted:
        t = Tree(treeFile)
    else:
        t = Tree(rootedMLtree)
    
    k = 1
    for n in t.traverse():
        if not n.is_root():
            if not n.is_leaf():
                n.add_features(name='n'+str(k))
                n.add_features(rootstrap=rootstrap_value[n.name]*100)
                k += 1
            else:
                n.add_features(rootstrap=rootstrap_value[n.name]*100)
    temp = os.path.splitext(treeFile)[0]+'.temp'
    rootstrapTree = os.path.splitext(treeFile)[0]+'.rootstrap'
    t.write(outfile=temp, features =["rootstrap"])
    x = dendropy.Tree.get(path=temp, schema='newick')
    x.write(path=rootstrapTree, schema='nexus')
    os.remove(temp)
    return polyphyly

def rootstrap():
    
    is_rooted = True
    out_group=None
    
    if len(sys.argv) == 1:
        raise SystemExit('Error: Please provide the ML tree file and the bootstrap trees file')
    elif len(sys.argv) == 2:
        if sys.argv[1] == '-h':
            print('\n')
            print('Syntax:\n')
            print('rootstrap.py <tree file> <bootstrap trees file> <is rooted> <outgroup file>\n')
            print('<treefile>  \t\t The tree file where you want to calculate the rootstrap support values in Newick format (e.g. tree.treefile)\n.')
            print('<bootstrap trees file>\t All the bootstrap trees in Newick format. e.g. .ufboot file from IQ-TREE (tree.ufboot)\n')
            print('<is rooted>  \t\t Are the ML tree and the boostrap trees rooted or not. default: True')
            print('\t\t\t True - The trees assumed to be rooted and no outgroup taxa infromation is required')
            print('\t\t\t False - The trees assumed to be unrooted and outgroup taxa infromation is required\n')
            print('<outgroup file> \t outgroup taxa in Nexus format.')
            print('\t\t\t The outgroup block can be part of the alignment file or in a separate file\n')
            return
        else:
            raise SystemExit('Error: unknown flag {}. Please provide the ML tree file and the bootstrap trees file'.format(sys.argv[1]))
    elif len(sys.argv) == 4:
        if sys.argv[3] == "False":
            raise SystemExit('Error: Please provide outgroup taxa in Nexus format')
    elif len(sys.argv) == 5:
        if sys.argv[3] == "False":
            is_rooted = False
            out_group = sys.argv[4]
    elif len(sys.argv) > 5:
        raise SystemExit('Error: too many arguments')
    print('Note: ingroup taxa are not monophyletic in {} bootstrap trees\n'.format(caluclate_rootstrap(sys.argv[1], sys.argv[2], is_rooted, out_group)))

if __name__ == '__main__':
    rootstrap()
