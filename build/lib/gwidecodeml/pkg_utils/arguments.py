#!/usr/bin/env python
# -*- coding: utf-8 -*-

from argparse import ArgumentParser
import os

args = None

def parse_arguments():
    parser = ArgumentParser(description="Parameters for GWideCodeml performance")
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument("-tree", dest="spp_tree", help="species tree in newick format", type=str)
    requiredNamed.add_argument("-work_dir", dest="wd", help="working_dir", type=str)
    requiredNamed.add_argument("-cds", dest="suffix", help="codon alignment files extension", type=str, default=".fa")
    parser.add_argument("-p", dest="threads", help="number of threads", type=str, default=1)
    parser.add_argument("-model", dest="mode", help="model testing", type=str, default="BS")
    parser.add_argument("-branch", dest="mark", help="text file containing species of the branch of interest", \
                        type=str, default=None)
    parser.add_argument("-no_lrt", dest="lrt_stop", action="store_true", default=False)
    parser.add_argument("-dnds", dest="get_dnds", action="store_true", default=False)
    parser.add_argument("-omin", dest="min_out", help="min nr. of outgroup spp", default=None)
    parser.add_argument("-cmin", dest="min_clade", help="min nr. of clade spp", default=None)
    parser.add_argument("-gene_trees", dest="single_trees",action="store_true", help="min nr. of outgroup spp", default=False)

    parser.parse_args()
    global args
    args = parser.parse_args()
    return args


args = parse_arguments()
global working_dir
working_dir = os.path.abspath(args.wd)

"""
# Trying required parameters
args = parse_arguments()
# define working dir
global working_dir
working_dir = os.path.abspath(args.wd)

# list codon alignments in working dir
suffix = args.suffix
fasta_files = [f for f in os.listdir(working_dir) if f.endswith(suffix)]
if len(fasta_files) == 0:
    print("Error, no files found with the provided extension. Exiting...")
    exit()

# Species tree
if args.spp_tree == None:
    print("Please, provide a tree in NEWICK format. Exiting...")
    exit()
else:
    spptree = os.path.abspath(args.spp_tree)

# Setting optional arguments:

# nr. threads define by the user
t = int(args.threads)

# Model selection
models = ["BS", "BM", "SM"]
model = str(args.mode)
if model not in models:
    print("Error, uncorrect model selection. Exiting...")
    exit()


elif model in ["BM", "BS"]:
    # Read branch marks
    if args.mark == None:
        print("If you choose branch or branch-site model, " \
              "a file specifying branch labels must be provided using -branch option. " \
              "Please, check the manual to provide a file in the right format.\n" \
              "Exiting... ")
        exit()
    else:
        branch_marks = read_branches(args.mark)

# branch labels
branch_file = args.mark

# lrt
lrt = args.lrt_stop

# dnds
dnds_file = args.get_dnds

# ougroup and clade min number
omin = args.min_out
cmin = args.min_clade
"""