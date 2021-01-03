#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import logging
import subprocess
import csv
import shutil
from Bio import SeqIO
from ete3 import EvolTree
from scipy import stats
from .arguments import args, working_dir


def list_files(work_dir, suffix):
    """Return a list of files with a certain pattern in a directory"""
    file_list = [f for f in os.listdir(work_dir) if f.endswith(suffix)]
    return file_list


def count_dict_values(list_in, dict_in, key):
    """Count nr values from a list having a given key in a dictionary"""
    total_values = 0
    for x in list_in:
        if x in dict_in[key]:
            total_values += 1
    return total_values


def fasta_ids(fasta):
    """Read multi-fasta file and return SPPTAG IDs in a list.
    Format >SPPTAG_geneid"""
    ids = []
    with open(fasta, "r", newline=None) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            ids.append(record.id.split("_")[0])
    return ids


def check_if_duplicates(list_elements):
    """Check if given list contains any duplicates"""
    if len(list_elements) == len(set(list_elements)):
        return False
    else:
        return True


def dup_tags(fasta_file_list):
    """Read a list of fasta files and returns fasta file name
    list containing duplicated genome tags"""
    fasta_dups = []
    tags = []
    for fasta_file in fasta_file_list:
        fasta_handle = SeqIO.parse(fasta_file, "fasta")
        tags = []
        for s in fasta_handle:
            genome_tag = s.id.split("_")[0]
            tags.append(genome_tag)
        # fasta_handle.close()
    if check_if_duplicates(tags):
        fasta_dups.append(fasta_file)
    return fasta_dups


def fasta2phy(msa_input, phy_out):
    """Convert fasta format into format compatible with paml"""
    input_handle = open(msa_input, "r", newline=None)
    output_handle = open(phy_out, "w")
    alignments = SeqIO.parse(input_handle, "fasta")
    msa_seqs = dict()
    for a in alignments:
        seq_id = str(a.id.split("_")[0])
        msa_seqs[seq_id] = str(a.seq)
        bases = len(str(a.seq))
    input_handle.close()
    # write first line of the output file
    output_handle.write("  " + str(len(msa_seqs.keys())) + "  " + str(bases) + "  " + "\n")
    # write sequences in paml format
    for k,v in msa_seqs.items():
        line_out = k + "  " + v + "\n"
        output_handle.write(line_out)
    output_handle.close()


def create_dir(work_dir, dir_name):
    """Create a new folder in the working_dir
    if folder exits, do nothing"""
    dirpath = os.path.join(work_dir, dir_name)
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        pass
        #shutil.rmtree(dirpath)
    else:
        os.mkdir(dirpath)


def move_files(work_dir,file_name,subfolder):
    """Move files to subfolder"""
    file_path = os.path.join(work_dir,file_name)
    shutil.move(file_path,subfolder)


def is_multiple_testing(branches_file):
    """Check if multiple foreground branches"""
    labels = branches_file.keys()
    multiple = 0
    for x in labels:
        if x not in ["0", "1"]:
            multiple += 1
    return multiple


def check_custom(files_path):
    """When custom model selected, checks if .ctl files are in dir
    and returns True if both in dir."""
    st1 = os.path.isfile(os.path.join(files_path, "custom_alt.ctl"))
    st2 = os.path.isfile(os.path.join(files_path, "custom_null.ctl"))
    if st1 and st2:
        return True
    else:
        return False


def prepare_codeml(round, fasta_file_name, species_tree, marks, wd):
    """Creates all files necessary for codeml performance,
    including marked, pruned tree, control files and alignment in phy format"""
    run_name = "r" + str(round).zfill(2)
    tree = EvolTree(species_tree)  # init tree every time a fasta is open
    name = fasta_file_name.replace(args.suffix, "")
    # print(name)
    # create path and change dir
    create_dir(wd, name)
    os.chdir(os.path.join(wd, name))
    genomes = fasta_ids(os.path.join(wd, fasta_file_name))  # genomes contained in fasta file
    # Tree prunning
    prune_tree(tree, genomes)

    # Mark branches if branch or branch-site models selected
    if args.mode in ["BM", "BS"]:
        mark_spp = list(set(marks[str(round)]).intersection(genomes))
        mark_branches(tree, mark_spp)

        # Check monophyly of taxa if -gene_trees option
        # Individual gene trees if -gene_trees option
        if args.single_trees:
            fast_tree(os.path.join(wd, fasta_file_name),
                      os.path.join(wd, name, fasta_file_name + ".ftree"))
            gene_tree = midpoint_root(tree_features(os.path.join(wd, name, fasta_file_name + ".ftree")))
            if not is_monophyletic(gene_tree, mark_spp):
                logging.warning("Check monophyly in the clade-of-interest: {}".format(name))

    tree.write(outfile=name + ".tree")  # write tree with only topology
    # File format converter: MSA fasta --> Phylip
    fasta2phy(os.path.join(wd, fasta_file_name), name + ".phy")

    # Create alt and null ctl files
    control_files(wd, args.mode, name, run_name)


def lnl(codeml_in):
    """Extract Ln Likelihood values from output files created by codeml"""
    with open(codeml_in, "r", newline=None) as file_in:
        if os.stat(codeml_in).st_size > 0:  # if file is not empty
            for line in file_in:
                if line.startswith("lnL"):
                    vals = re.findall(r"[-]?\d*\.\d+|\d+", line)
                    for x in vals:
                        if float(x) < 0:
                            log = float(x)
                    break
        else:
            log = None
    return log


def chisqprob(chisq, df):
    """Calculate p-vale of chis-square"""
    value = stats.chi2.sf(chisq, df)
    return value


def lrt(ln_1, ln_2, df=1):
    """Likelihood Ratio Test calculator
    degrees of freedom (df, default = 1)"""
    if ln_1 and ln_2:
        val = 2 * (float(ln_1) - (float(ln_2)))
        p_val = chisqprob(val, df)
    else:
        p_val = None
    return p_val


def model_selection(select):
    """Create dicts with model parameters"""
    # branch-site model null hypothesis
    bs0 = {
        "model": "2",
        "NSsites": "2",
        "fix_omega": "1",
        "omega": "1",
        "ncatG": "3"
    }
    # branch-site model alt hypothesis
    bs1 = {
        "model": "2",
        "NSsites": "2",
        "fix_omega": "0",
        "omega": "2",
        "ncatG": "3"
    }
    # branch model null hypothesis
    bm0 = {
        "model": "0",
        "NSsites": "0",
        "fix_omega": "0",
        "omega": "2",
        "ncatG": "3"
    }
    # branch model alt hypothesis
    bm1 = {
        "model": "2",
        "NSsites": "0",
        "fix_omega": "0",
        "omega": "2",
        "ncatG": "3"
    }
    # site model M1a null hypothesis
    sm0 = {
        "model": "0",
        "NSsites": "1",
        "fix_omega": "0",
        "omega": "2",
        "ncatG": "3"
    }
    # site model M2a alt hypothesis
    sm1 = {
        "model": "0",
        "NSsites": "2",
        "fix_omega": "0",
        "omega": "2",
        "ncatG": "3"
    }
    if select == "BS":
        return bs0, bs1
    elif select == "BM":
        return bm0, bm1
    elif select == "SM":
        return sm0, sm1


def codeml_settings(seq_input, tree_input, out_name, model):
    """Create and edit control files to run codeml according to model selected"""

    ctl_find = ["seq_file.phy", "tree.nwk", "out_name", "model = 2", "NSsites = 2",
                "fix_omega = 0", "omega = .4", "ncatG = 3"]
    ctl_replace = [seq_input, tree_input, out_name]
    ctl_replace.append("model = " + model["model"])
    ctl_replace.append("NSsites = " + model["NSsites"])
    ctl_replace.append("fix_omega = " + model["fix_omega"])
    ctl_replace.append("omega = " + model["omega"])
    ctl_replace.append("ncatG = " + model["ncatG"])

    return ctl_find, ctl_replace


def custom_settings(seq_input, tree_input, out_name):
    """Create and edit control files for custom models"""
    ctl_find = ["seq_file.phy", "tree.nwk", "out_name"]
    ctl_replace = [seq_input, tree_input, out_name]

    return ctl_find, ctl_replace


def read_replace(f_in, f_out, findlines, replacelines):
    """Replace multiple lines in a text file"""
    find_replace = dict(zip(findlines, replacelines))
    try:
        with open(f_in, "r", newline=None) as data:
            with open(f_out, 'w') as new_data:
                for line in data:
                    for key in find_replace:
                        if key in line:
                            line = line.replace(key, find_replace[key])
                    new_data.write(line)
    except FileNotFoundError:
        logging.error("Error: gwcodeml.ctl. No such file in the working directory")
        exit()


def control_files(work_dir, model, gene, run):
    """Creates both null and alt ctl files"""
    if model != "custom":
        # Create control (.ctl) files
        ctl_in = os.path.join(work_dir, "gwcodeml.ctl")
        h0, h1 = model_selection(model)
        # create ctl file for null hypothesis testing
        f0, r0 = codeml_settings(gene + ".phy", gene + ".tree", gene + "_" + run + "_null.txt", h0)
        read_replace(ctl_in, gene+"_null.ctl", f0, r0)
        # create ctl file for alternative hypothesis testing
        f1, r1 = codeml_settings(gene + ".phy", gene + ".tree", gene + "_" + run + "_alt.txt", h1)
        read_replace(ctl_in, gene + "_alt.ctl", f1, r1)
    elif model == "custom":
        # create _null.ctl
        ctl_null = os.path.join(work_dir, "custom_null.ctl")
        f0, r0 = custom_settings(gene + ".phy", gene + ".tree", gene + "_" + run + "_null.txt")
        read_replace(ctl_null, gene + "_null.ctl", f0, r0)
        # create _alt.ctl
        ctl_alt = os.path.join(work_dir, "custom_alt.ctl")
        f1, r1 = custom_settings(gene + ".phy", gene + ".tree", gene + "_" + run + "_alt.txt")
        read_replace(ctl_alt, gene + "_alt.ctl", f1, r1)


def run_codeml(ctlfile):
    """Create individual folders for running codeml"""
    path_to_gene = os.path.join(working_dir, ctlfile)
    os.chdir(path_to_gene)
    # codeml null hypothesis
    logging.info("Running codeml null hypothesis on {}".format(ctlfile))
    bash_command = "codeml " + ctlfile + "_null.ctl"
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    # codeml alt hypothesis
    logging.info("Running codeml alternative hypothesis on {}".format(ctlfile))
    bash_command = "codeml " + ctlfile + "_alt.ctl"
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    return output, error


def omega_bs(alt_file):
    """Extract omega values from branch-site results"""
    with open(alt_file, "r", newline=None) as alt:
        lines = alt.readlines()
    for line in lines:
        if line.startswith("foreground w"):
            found = line.replace("\n", "")
            omegas = re.findall("\d+\.\d+", found)
            return float(omegas[-1])


def omega_bm(alt_file):
    """Extract omega values from branch results"""
    with open(alt_file, "r", newline=None) as alt:
        lines = alt.readlines()
    for line in lines:
        if line.startswith("w (dN/dS) for branches:"):
            found = line.replace("\n", "")
            omegas = re.findall("\d+\.\d+", found)
            return omegas


def beb(alt_file):
    """Extract significant aa positions according to BEB test
    for site and branch-site models"""
    flist = open(alt_file).readlines()

    aas = dict()
    parsing = False
    for line in flist:
        if line.startswith("Bayes Empirical Bayes"):
            parsing = True
        elif line.startswith("The grid "):
            parsing = False
        if parsing:
            text = line.strip()
            if re.match(r"^\d+.*$",text):
                text_list = [x for x in text.split(" ") if x != ""]
                if "*" in text_list[2]:
                    aas[text_list[0]] = text_list[1]
    return aas


def output_bs(gene_name, alt_file):
    """Return significant aa positions and amino-acid changes for branch-site analysis"""
    beb_dict = beb(alt_file)
    for k, v in beb_dict.items():
        if v != None:
            l_out = gene_name + "\t" + k + "\t" + v
            return l_out


def output_sm(gene_name, alt_file):
    """Return significant aa positions and amino-acid changes for site analysis"""
    beb_dict = beb(alt_file)
    nr_pos = 0
    pos = []
    for k,v in beb_dict.items():
        if v != None:
            nr_pos += 1
            pos.append(k)
    if len(pos) != 0:
        line_out = gene_name + "\t" + str(nr_pos)+ "\t" + ";".join(pos)
        return line_out


def output_bm(gene_name, alt_file):
    """Return significant branches with omegas > 1 in branch analysis"""
    omegas = omega_bm(alt_file)
    for c, value in enumerate(omegas[1:],1):
        if float(value) > 1:
            line_out = gene_name +"\t"+"#"+str(c)
            return line_out


def if_significant(alt_file, null_file, model):
    """Returns true when a gene analyzed
    has a pval < 0.05 according to LRT"""
    ln_alt = lnl(alt_file)
    ln_null = lnl(null_file)
    p_val = None
    if model in ["BS", "BM"]:
        p_val = lrt(ln_alt, ln_null)
    elif model == "SM":
        # in site-model, LRT has 2 degrees of freedom
        p_val = lrt(ln_alt, ln_null, 2)
    if p_val < 0.05:
        return True
    else:
        return False


def final_output(sign_list, model, run, dirpath):
    """Write significant genes to a tsv output file"""
    out_file_name = "results_" + run + "_" + model + ".tsv"
    out_file = open(os.path.join(dirpath, out_file_name), "w")

    # branch model: get omega and end
    if model == "BM":
        out_file.write("Gene_name" + "\t" + "Branch_label" + "\n")
        for gene in sign_list:
            alt_file = os.path.join(dirpath, gene, gene + "_" + run + "_alt.txt")
            l_out = output_bm(gene, alt_file)
            if l_out:
                out_file.write(l_out + "\n")

    # branch-site model: omega + beb positions
    elif model == "BS":
        out_file.write("Gene_name" + "\t" + "AA_position" + "\t" + "AA" + "\n")
        for gene in sign_list:
            alt_file = os.path.join(dirpath, gene, gene+"_"+run+"_alt.txt")
            o = omega_bs(alt_file)
            if float(o) > 1:
                l_out = output_bs(gene, alt_file)
                if l_out:
                    out_file.write(l_out+"\n")

    # site model: beb positions under PS
    elif model == "SM":
        out_file.write("Gene_name" + "\t" + "Nr_positions" + "\t" + "AA Positions" + "\n")
        for gene in sign_list:
            alt_file = os.path.join(dirpath, gene, gene + "_" + run + "_alt.txt")
            l_out = output_sm(gene, alt_file)
            if l_out:
                out_file.write(l_out+"\n")

    out_file.close()


def dnds_output(model, work_dir, genes_list, run):
    """Creates an extra file extracting all dnds values from all analyzed genes"""
    omega_file_name = "dnds_" + run + "_" + model + ".tsv"
    omega_file = open(omega_file_name, "w")
    tsvwriter = csv.writer(omega_file, delimiter="\t")
    if model == "BM":
        col_names = ["Gene","#0","#1"]
        tsvwriter.writerow(col_names)
        for gene in genes_list:
            dnds_values = omega_bm(os.path.join(work_dir, gene, gene + "_" + run + "_alt.txt"))
            l_out = [gene] + dnds_values
            tsvwriter.writerow(l_out)
    elif model == "BS":
        col_names = ["Gene", "#1"]
        tsvwriter.writerow(col_names)
        for gene in genes_list:
            dnds_branch = omega_bs(os.path.join(work_dir, gene, gene + "_" + run + "_alt.txt"))
            l_out = [gene, dnds_branch]
            tsvwriter.writerow(l_out)
    omega_file.close()


def read_branches(txt_in):
    """Read branch labels from a text file
    return info as dict"""
    branches = dict()
    with open(txt_in, "r", newline=None) as f_in:
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


def node_dict(tree_in):
    """Parse tree and return a dictionary with leaf name as keys
    and node_ids as values"""
    nodes = dict()
    for leaf in tree_in:
        nodes[leaf.name] = str(leaf.node_id)
    return nodes


def prune_tree(ptree,labels):
    """Prune tree until it contains labels provided"""
    nodes = read_leaves(ptree)
    ptree.prune(list(set(nodes).intersection(labels)))


def mark_branches(mtree, labels):
    tree_nodes = node_dict(mtree)
    if len(labels) > 1:
        ancestor = mtree.get_common_ancestor(labels)
        nx = str(ancestor.node_id)
        mtree.mark_tree([nx], marks=["#1"])
    elif len(labels) == 1:
        n1 = tree_nodes[labels[0]]
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

