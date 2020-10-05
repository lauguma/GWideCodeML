#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
from ete3 import EvolTree
import subprocess

def read_branches(txt_in):
    """Read branch labels from a text file
    return info as dict"""
    branches = dict()
    with open(txt_in, "r",newline=None) as f_in:
        freader = csv.reader(f_in, delimiter=",")
        for row in freader:
            if row[1] in branches:
                branches[row[1]].append(row[0])
            else:
                branches[row[1]] = [row[0]]

    return branches

def read_leaves(tree_in):
    """Parse species tree and return a list with leaves names"""
    # save in a list all spp contained in spp tree
    leaves = []
    for leaf in tree_in:
        leaves.append(leaf.name)
    return leaves

def prune_tree(ptree,labels):
    """Prune tree until it contains labels provided"""
    nodes = read_leaves(ptree)
    ptree.prune(list(set(nodes).intersection(labels)))

def mark_branches(mtree, labels):
    if len(labels) > 1:
        ancestor = mtree.get_common_ancestor(labels)
        nx = str(ancestor.node_id)
        mtree.mark_tree([nx], marks=["#1"])
    else:
        n1 = str(labels[0].node_id)
        mtree.mark_tree([n1], marks=["#1"])



def is_monophyletic(tree_object, tag_names):
    """Check whether a group of individuals group exclusively
    together within a tree partition"""
    t = tree_object.check_monophyly(values=tag_names, target_attr="spptag")[0]
    if t == True:
        return True
    else:
        return False

def fast_tree(ali_file,ali_out):
    """Run FastTree to get a phylogeny per gene"""
    ft_command = "FastTree -gtr -nt " + ali_file
    process = subprocess.Popen(ft_command.split(), stdout=subprocess.PIPE)
    process_output,error = process.communicate()
    file = open(ali_out, "w")
    file.write(process_output.decode('utf-8'))
    file.close()

def tree_features(tree_file):
    """Add a feature to the tree called <spptag> with the Species Tag"""
    tree_handle = EvolTree(tree_file)
    for leaf in tree_handle.iter_leaves():
        leaf.add_feature("spptag", leaf.name.split("_")[0])
    return tree_handle

def midpoint_root(genetree):
    """Gene trees will be temporary midpoint rooted to check monophyly"""
    r = genetree.get_midpoint_outgroup()
    genetree.set_outgroup(r)
    return genetree

