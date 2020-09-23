# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 14:12:29 2019

@author: Suha Naser-Khdour
"""

import sys
from ete3 import Tree

def tree_with_internal_names(treeFile):
    t = Tree(treeFile)
    k = 1
    for n in t.traverse():
        if not n.is_root():
            if not n.is_leaf():
                n.add_features(name='n'+str(k))
                k += 1
        else:
            n.add_features(name='n0')
    print(t.get_ascii(show_internal=True))
    return 

def dist_from_true(treefile):
    tree_with_internal_names(treefile)
    trueRoot = input('where do you believe the true root to be? (choose the node at the end of that branch)\n')
    t = Tree(treefile)
    k = 1
    for n in t.traverse():
        if not n.is_root():
            if not n.is_leaf():
                n.add_features(name='n'+str(k))
                k += 1
        else:
            n.add_features(name='n0')
    if trueRoot != 'n0':
        trueRoot = t&trueRoot
        distFromRoot = t.get_distance(trueRoot, topology_only=True)
        startPonit = t.get_distance(trueRoot)
        endPoint = t.get_distance(trueRoot.up)
        print('rSED = {}\n'.format(distFromRoot))
        print('rBED = {} - {}\n'.format(endPoint,startPonit))
    else:
        trueRoot = t&trueRoot
        distFromRoot = t.get_distance(trueRoot, topology_only=True)
        startPonit = 0.0
        endPoint = 0.0
        for n in trueRoot.children:
            if t.get_distance(n) > endPoint:
                endPoint = t.get_distance(n)
        print('rSED = {}\n'.format(distFromRoot))
        print('rBED = 0.0 - {}\n'.format(endPoint))
        
    #option 1: the distance from the midpoint of the ML root branch to the true root branch
    branchMidpoint = t.get_distance(trueRoot, trueRoot.up)/2
    distFromRoot = t.get_distance(trueRoot) - branchMidpoint
    #option 2: the number of branches between the ML root branch and the true root branch
    return [distFromRoot,startPonit,endPoint]

def rBED_rSED():
    if len(sys.argv) == 1:
        raise SystemExit('Error: Please provide the rooted ML tree file')
    elif len(sys.argv) > 2:
        raise SystemExit('Error: too many arguments')
    dist_from_true(sys.argv[1])

if __name__ == '__main__':
    #required input
    rBED_rSED()            
