#!/usr/bin/env python
# -*- coding: utf-8 -*-

from argparse import ArgumentParser
import os


def parse_arguments():
    parser = ArgumentParser(description="Parameters for GWideCodeml performance", prog="GWideCodeML")
    required_named = parser.add_argument_group('required arguments')
    required_named.add_argument("-tree", dest="spp_tree", help="species tree in Newick format", type=str)
    required_named.add_argument("-work_dir", dest="wd", help="working_dir", type=str)
    required_named.add_argument("-cds", dest="suffix", help="codon alignment files extension", type=str, default=".fa")
    parser.add_argument("-p", dest="threads", help="maximum number of threads", type=str, default=1)
    parser.add_argument("-model", dest="mode", help="model testing", type=str, default="BS")
    parser.add_argument("-branch", dest="mark", help="text file containing species labels: "
                                                     "clade-of-interest (1,2...) vs outgroups (0)",
                        type=str, default=None)
    parser.add_argument("-no_lrt", dest="lrt_stop", help="do not perform LRTs", action="store_true", default=False)
    parser.add_argument("-dnds", dest="get_dnds",
                        help="obtain an additional file with dN/dS values for all genes (only BM/BS models)",
                        action="store_true", default=False)
    parser.add_argument("-omin", dest="min_out", help="min nr. of outgroup spp", default=None)
    parser.add_argument("-cmin", dest="min_clade", help="min nr. of clade-of-interest spp", default=None)
    parser.add_argument("-gene_trees", dest="single_trees", action="store_true",
                        help="individual gene trees are created with FastTree and used for check monophyly "
                             "in the foreground branch", default=False)
    parser.add_argument('--version', action='version', version='%(prog)s 1.1', help="prints current version")

    parser.parse_args()
    global args
    args = parser.parse_args()
    return args


args = parse_arguments()
global working_dir
working_dir = os.path.abspath(args.wd)
