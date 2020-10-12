#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import multiprocessing as mp
from ete3 import EvolTree
from .pkg_utils import codemltools as ct
from .pkg_utils import treetools as tt
from .pkg_utils.arguments import args
from .pkg_utils import utils
import logging


# Main steps of the program
def main():

    working_dir = os.path.abspath(args.wd)
    # logging
    log_file = os.path.join(working_dir, "gwidecodeml.log")
    logging.basicConfig(filename=log_file, filemode="w", level=logging.INFO,
                        format='%(asctime)s :: %(levelname)s :: %(message)s')
    logging.info("Working directory: {}".format(working_dir))

    # list codon alignments in working dir
    fasta_files = utils.list_files(working_dir, args.suffix)
    if len(fasta_files) == 0:
        logging.error("No alignments found under the provided extension. Exiting...")
        exit()

    # Species tree
    if args.spp_tree == None:
        logging.error("Please, provide a tree in NEWICK format. Exiting...")
        exit()
    spptree = os.path.abspath(args.spp_tree)


    # Setting optional arguments:
    # nr. threads define by the user
    t = int(args.threads)

    # Model selection
    models = ["BS", "BM", "SM", "custom"]
    model = str(args.mode)
    if model not in models:
        logging.error("Incorrect model selection. Exiting...")
        exit()
    elif model == "custom":
        if not utils.check_custom(working_dir):
            logging.error("Custom control files for codeml not found. Please, create and copy"
                          "custom_alt.ctl and custom_null.ctl files to the working directory to "
                          "run GWideCodeML under a custom mode. Exiting...")
            exit()

    rounds = 1  # one round by default unless multiple testing
    if model in ["BM", "BS"]:
        if args.mark == None:
            logging.error("If you choose branch or branch-site model, a file specifying \
            branch labels must be provided using -branch option.\
            Please, check the manual to provide a file in the right format.\n Exiting... ")
            exit()
        else:
            branch_marks = tt.read_branches(args.mark)
            # check if multiple branch testing
            if utils.is_multiple_testing(branch_marks) != 0:
                rounds += utils.is_multiple_testing(branch_marks)
                logging.info("Multiple branch testing selected: {} branches will be "
                             "analyzed as foreground branches.".format(str(rounds)))

    logging.info("Total number of alignments provided: {}".format(str(len(fasta_files))))
    logging.info("Model selected: {}".format(model))
    logging.info("Species tree: {}".format(spptree))
    logging.info("Nr. of threads: {}".format(str(t)))

    # START
    outputs = []
    for r in range(1, rounds+1):
        run_name = "r" + str(r).zfill(2)
        logging.info("Starting GWideCodeML analysis... ROUND NAME: {}".format(run_name))
        # Filter fasta if contains duplicated genome tags
        remove_fasta = []
        fasta_dups = utils.dup_tags(fasta_files)
        remove_fasta.extend(fasta_dups)
        if len(fasta_dups) > 0:
            logging.warning("Found duplicated genome tags in fasta files. {} fasta files "
                            "won't be analyzed".format(len(fasta_dups)))
            logging.warning("Discarded fasta files: {}".format(fasta_dups))
        # Filter fasta files if omin and cmin provided
        if args.min_out or args.min_clade:
            logging.info("Applying filters...")
            # If omin or cmin, branch marks file is necessary
            if args.mark == None:
                logging.error("If you choose filter out alignments by a minimum number of taxa \
                and/or outgroups, you should provide labels in a text file using the -branch option. \
                Please, check the manual to provide a file in the right format.\n \
                Exiting... ")
                exit()
            else:
                branch_marks = tt.read_branches(args.mark)

                if args.min_out and branch_marks:
                    logging.info("Minimum number of outgroup species: {}".format(str(args.min_out)))
                    c = 0
                    for s in fasta_files:
                        spp_tags = utils.fasta_ids(os.path.join(working_dir, s))
                        outgroups = utils.count_dict_values(spp_tags, branch_marks, "0")
                        if outgroups < int(args.min_out):
                            remove_fasta.append(s)
                            c += 1
                    logging.info("{} alignments removed by min. outgroup filter.".format(str(c)))

                if args.min_clade and branch_marks:
                    logging.info("Minimum number of clade-of-interest species: {}".format(str(args.min_clade)))
                    c = 0
                    for s in fasta_files:
                        spp_tags = utils.fasta_ids(os.path.join(working_dir, s))
                        clade = utils.count_dict_values(spp_tags, branch_marks, str(r))
                        if clade < int(args.min_clade):
                            remove_fasta.append(s)
                            c += 1
                    logging.info("{} alignments removed by min. clade-of-interest filter.".format(str(c)))

        if len(remove_fasta) > 0:
            logging.info("These alignments won't be analyzed: ")
        for item in remove_fasta:
            logging.info(item)
        # Filter fastas if necessary
        fasta_list = [x for x in fasta_files if x not in remove_fasta]

        logging.info("Starting analysis of {} alignments".format(str(len(fasta_list))))
        # Prepare files that have passed filters for codeml performance
        gene_names = []
        for fasta in fasta_list:
            tree = EvolTree(spptree)  # init tree every time a fasta is open
            name = fasta.replace(args.suffix,"")
            #print(name)
            # create path and change dir
            utils.create_dir(working_dir, name)
            os.chdir(os.path.join(working_dir, name))
            genomes = utils.fasta_ids(os.path.join(working_dir, fasta)) # genomes contained in fasta file
            gene_names.append(name)
            # Tree prunning
            tt.prune_tree(tree,genomes)
            # Individual gene trees
            tt.fast_tree(os.path.join(working_dir, fasta), os.path.join(working_dir, name, fasta + ".ftree"))
            gene_tree = tt.midpoint_root(tt.tree_features(os.path.join(working_dir, name, fasta + ".ftree")))

            # Mark branches if branch or branch-site models selected
            if model in ["BM", "BS"]:
                mark_spp = list(set(branch_marks[str(r)]).intersection(genomes))
                tt.mark_branches(tree,mark_spp)
                # logging.info("Foreground branch is composed of {} labels".format(mark_spp))
                # Check monophyly of taxa
                if not tt.is_monophyletic(gene_tree, mark_spp):
                    logging.warning("Check monophyly in the clade-of-interest: {}".format(name))

            tree.write(outfile=name+".tree") # write tree with only topology
            # File format converter: MSA fasta --> Phylip
            utils.fasta2phy(os.path.join(working_dir, fasta), name + ".phy")

            # Create alt and null ctl files
            ct.control_files(working_dir, args.mode, name, run_name)

        logging.info("Control files successfully created, GWidecodeml is ready for codeml performance")

        # Run codeml in parallel
        pool = mp.Pool(t) # number of threads
        pool.map(ct.run_codeml, [ctl for ctl in gene_names])
        pool.close()
        os.chdir(working_dir)

        if args.lrt_stop:
            logging.info("Codeml has finished, output files successfully created. Exiting...")
            exit()
        else:
            logging.info("Codeml has finished, output files successfully created. Running LRTs...")


        # LRTs, keep genes with p-value < 0.05
        significants = []
        for x in gene_names:
            hyp_alt = os.path.join(working_dir, x, x + "_" + run_name + "_alt.txt")
            hyp_null = os.path.join(working_dir, x, x + "_" + run_name + "_null.txt")
            if ct.if_significant(hyp_alt, hyp_null, args.mode):
                significants.append(x)

        # Write significant genes to output
        ct.final_output(significants, args.mode, run_name, working_dir)
        logging.info("Total nr. of genes rejecting null hypothesis: {}".format(len(significants)))

        logging.info("dnds option selected. Omega values will be written to an output file.")
        # if -dnds option, an extra file is created with omega values for all tested genes
        if args.get_dnds and args.mode in ["BM","BS"]:
            ct.dnds_output(args.mode, working_dir, gene_names, run_name)



        out_file_name = "results_{}_{}.tsv".format(run_name,args.mode)
        outputs.append(out_file_name)
        logging.info("Round name {} finished".format(run_name))

    logging.info("GWideCodeML successfully run. Please, check results files: {}".format(outputs))



if __name__ == "__main__":
    main()

