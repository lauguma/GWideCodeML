#!/usr/bin/python

'''
GWideCodeML: a Python automatic pipeline for running codeml
'''

__author__ = "Laura Gutierrez Macias and Christina Toft"
__version__ = "1.0.1"
__mantainer__ = "Laura G. Macias"
__email__ = "laugmacias@gmail.com"


import os
import subprocess
from argparse import ArgumentParser
import multiprocessing as mp
import csv
from scipy import stats
from ete3 import EvolTree
from Bio import SeqIO
import re


# Read branches of interest from a text file
# Return as a dictionary
def read_branches(txt_in):
    branches = dict()
    with open(txt_in, "r") as f_in:
        freader = csv.reader(f_in,delimiter = ",")
        for row in freader:
            if row[1] in branches:
                branches[row[1]].append(row[0])
            else:
                branches[row[1]] = [row[0]]

    return branches

# Count nr values from a list having a given key in a dictionary
def count_dict_values(list_in,dict_in,key):
    total_values = 0
    for x in list_in:
        if x in dict_in[key]:
            total_values += 1
    return total_values

# Parse species tree and return a list with leaves names
def read_leaves(tree_file):
    tree_handle = EvolTree(tree_file)
    leaves = [] # save in a list all spp contained in spp tree
    for leaf in tree_handle:
        leaves.append(leaf.name)
    return leaves

# Read multi-fasta file and return sequence IDs in a list
# Format spptag_geneid. Get only spp tag
def fasta_ids(fasta):
    ids = []
    fasta_file = SeqIO.parse(fasta, "fasta")
    for s in fasta_file:
        ids.append(s.id.split("_")[0])
    fasta_file.close()
    return ids

# Transform fasta format to format compatible with paml
def fasta2phy(msa_input, phy_out):
    input_handle = open(msa_input, "rU")
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

# Extract Ln Likelihood values from output files created by codeml
def lnl(codeml_in):
    with open(codeml_in, "r") as file_in:
        if os.stat(codeml_in).st_size > 0:  # if file is not empty
            for line in file_in:
                if line.startswith("lnL"):
                    vals = re.findall(r"[-]?\d*\.\d+|\d+", line)
                    for x in vals:
                        if float(x) < 0:
                            lnl = float(x)
                    break

        else:
            lnl = None
    return lnl

# Likelihood ratio test with likelihood values from null and alternative hypotheses
# df: degrees of freedom, default=1. In site model is df=2
def lrt(ln_1, ln_2,df=1):
    stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)
    if ln_1 and ln_2:
        val = 2 * (float(ln_1) - (float(ln_2)))
        p_val = stats.chisqprob(val, df)
    else:
        p_val = None

    return p_val

# Create dictionaries containing model parameters
def model_selection(select):
    # branch-site model null hypothesis
    bs0 = {
        "model":"2",
        "NSsites":"2",
        "fix_omega":"1",
        "omega":"1",
        "ncatG":"3"
    }
    # branch-site model alt hypothesis
    bs1 = {
        "model":"2",
        "NSsites":"2",
        "fix_omega":"0",
        "omega":"2",
        "ncatG":"3"
    }
    # branch model null hypothesis
    bm0 = {
        "model":"0",
        "NSsites":"0",
        "fix_omega":"0",
        "omega":"2",
        "ncatG":"3"
    }
    # branch model alt hypothesis
    bm1 = {
        "model":"2",
        "NSsites":"0",
        "fix_omega":"0",
        "omega":"2",
        "ncatG":"3"
    }
    # site model M1a null hypothesis
    sm0 = {
        "model":"0",
        "NSsites":"1",
        "fix_omega":"0",
        "omega":"2",
        "ncatG":"3"
    }
    # site model M2a alt hypothesis
    sm1 = {
        "model":"0",
        "NSsites":"2",
        "fix_omega":"0",
        "omega":"2",
        "ncatG":"3"
    }

    if select == "BS":
        return bs0,bs1
    elif select == "BM":
        return bm0,bm1
    elif select == "SM":
        return sm0,sm1

# replace multiple lines in a text file
def read_replace(f_in,f_out,findlines,replacelines):
    find_replace = dict(zip(findlines, replacelines))
    try:
        with open(f_in,"r") as data:
            with open(f_out, 'w') as new_data:
                for line in data:
                    for key in find_replace:
                        if key in line:
                            line = line.replace(key, find_replace[key])
                    new_data.write(line)
    except FileNotFoundError:
        print("Error: gwcodeml.ctl. No such file in the working directory")
        exit()


# Create and edit control files to run codeml according to model selected
def codeml_settings(seq_input, tree_input, out_name, model):

    ctl_find = ["seq_file.phy","tree.nwk","out_name","model = 2",\
                "NSsites = 2","fix_omega = 0","omega = .4","ncatG = 3"]
    ctl_replace = [seq_input,tree_input,out_name]

    ctl_replace.append("model = " + model["model"])
    ctl_replace.append("NSsites = " + model["NSsites"])
    ctl_replace.append("fix_omega = " + model["fix_omega"])
    ctl_replace.append("omega = " + model["omega"])
    ctl_replace.append("ncatG = " + model["ncatG"])

    return ctl_find,ctl_replace


# For codeml running, create individual folders
# located at the working dir to run codeml on them, both null and alt
def run_codeml(ctlFile):
    os.chdir(os.path.join(working_dir,ctlFile))
    # codeml null hypothesis
    bash_command = "codeml " + ctlFile + "_null.ctl"
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    # codeml alt hypothesis
    bash_command = "codeml " + ctlFile + "_alt.ctl"
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    return output, error

# omega for branch-site models
def omega_bs(alt_file):
    with open(alt_file, "r") as alt:
        lines = alt.readlines()
    for line in lines:
        if line.startswith("foreground w"):
            found = line.replace("\n", "")
            omegas = re.findall("\d+\.\d+", found)
            return float(omegas[-1])

# omega for branch model (one omega for the whole coding seq)
def omega_bm(alt_file):
    with open(alt_file, "r") as alt:
        lines = alt.readlines()
    for line in lines:
        if line.startswith("w (dN/dS) for branches:"):
            found = line.replace("\n", "")
            omegas = re.findall("\d+\.\d+", found)

            return omegas


# get beb positions for site models, return a dictionary
def beb(alt_file):
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

def output_bs(gene_name):
    beb_dict = beb(gene_name + "_alt.txt")
    nr_pos = 0
    pos = []
    for k,v in beb_dict.items():
        if v != None:
            l_out = gene_name + "\t" + k + "\t" + v
            return l_out

def output_sm(gene_name):
    beb_dict = beb(gene_name + "_alt.txt")
    nr_pos = 0
    pos = []
    for k,v in beb_dict.items():
        if v != None:
            nr_pos += 1
            pos.append(k)
    if len(pos) != 0:
        line_out = gene_name + "\t" + str(nr_pos)+ "\t" + ";".join(pos)
        return line_out

def output_bm(gene_name):
    omegas = omega_bm(gene_name + "_alt.txt")
    for c, value in enumerate(omegas[1:],1):
        if float(value) > 1:
            line_out = gene_name +"\t"+"#"+str(c)
            return line_out


# Main steps of the program
def main():
    parser = ArgumentParser(description="Parameters for GWideCodeml performance")
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument("-tree", dest="spp_tree", help="species tree in newick format", type=str)
    requiredNamed.add_argument("-work_dir", dest="wd", help="working_dir", type=str)
    requiredNamed.add_argument("-cds", dest="suffix", help="codon alignment files extension", type=str, default=".fa")
    parser.add_argument("-p", dest="threads", help="number of threads", type=str,default=1)
    parser.add_argument("-model", dest="mode", help="model testing", type=str, default="BS")
    parser.add_argument("-branch", dest="mark", help="text file containing species of the branch of interest", type=str)
    parser.add_argument("-no_lrt", dest="lrt_stop", action="store_true", default=False)
    parser.add_argument("-omin", dest="min_out", help="min nr. of outgroup spp", default=False)
    parser.add_argument("-cmin", dest="min_clade", help="min nr. of clade spp", default=False)

    args = parser.parse_args()

    # Trying required parameters

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
    spptree = os.path.abspath(args.spp_tree)

    # Setting optional arguments:

    # nr. threads define by the user
    t = int(args.threads)

    # Model selection
    models = ["BS","BM","SM"]
    model = str(args.mode)
    if model not in models:
        print("Error, uncorrect model selection. Exiting...")
        exit()


    elif model in ["BM","BS"]:
        # Read branch marks
        if args.mark == None:
            print("If you choose branch or branch-site model, " \
                  "a file specifying branch labels must be provided using -branch option. " \
                  "Please, check the manual to provide a file in the right format.\n" \
                  "Exiting... ")
            exit()
        else:
            branch_marks = read_branches(args.mark)

    
    # Filter fasta files if omin and cmin provided
    remove_fasta = []
    if args.min_out:
        for s in fasta_files:
            spp_tags = fasta_ids(os.path.join(working_dir, s))
            outgroups = count_dict_values(spp_tags, branch_marks, "0")
            if outgroups < int(args.min_out):
                remove_fasta.append(s)

    if args.min_clade:
        for s in fasta_files:
            spp_tags = fasta_ids(os.path.join(working_dir, s))
            clade = count_dict_values(spp_tags, branch_marks, "1")
            if clade < int(args.min_clade):
                remove_fasta.append(s)



    # Filter fastas if necessary
    fasta_list = [x for x in fasta_files if x not in remove_fasta]

    # Prepare files that have passed filters for codeml performance
    gene_names = []
    for fasta in fasta_list:
        tree = EvolTree(spptree)  # init tree every time a fasta is open
        nodes = read_leaves(spptree)
        name = fasta.replace(suffix,"")
        #print(name)
        os.mkdir(os.path.join(working_dir,name)) # create one folder per gene analyzed
        os.chdir(os.path.join(working_dir,name))
        genomes = fasta_ids(os.path.join(working_dir, fasta)) # genomes contained in fasta file
        gene_names.append(name)

        # Tree prunning
        spp2tree = list(set(nodes).intersection(genomes)) # spp to keep in tree
        tree.prune(spp2tree)
        # Mark branches if branch or branch-site models selected
        # branch site model, one branch each time
        if model == "BS":
            mark_spp = branch_marks["1"]
            ancestor = tree.get_common_ancestor(mark_spp)  # get ancestor of spp to mark branches
            n1 = str(ancestor.node_id)
            tree.mark_tree([n1], marks=["#1"])
            tree.write(outfile=name+".tree") # write tree with only topology

        # branch model allows for multiple branch marks
        elif model == "BM":
            nodes = []
            hashes = []
            for x in branch_marks.keys():
                if x != "0": # if not outgroup
                    ancestor = tree.get_common_ancestor(branch_marks[x])  # get ancestor of spp to mark branches
                    n1 = str(ancestor.node_id)
                    nodes.append(n1)
                    hashes.append("#"+x)
            tree.mark_tree(nodes, marks=hashes)

        tree.write(outfile=name+".tree") # write tree with only topology

        # File format converter: MSA fasta --> Phylip
        fasta2phy(os.path.join(working_dir,fasta), name+".phy")


    # step 4: create .ctl files
        ctl_in = os.path.join(working_dir,"gwcodeml.ctl")
        h0,h1 = model_selection(model)
        # create ctl file for null hypothesis testing
        f0,r0 = codeml_settings(name+".phy", name+".tree", name+"_null.txt", h0)
        read_replace(ctl_in, name+"_null.ctl", f0, r0)

        # create ctl file for alternative hypothesis testing
        f1, r1 = codeml_settings(name + ".phy", name + ".tree", name + "_alt.txt", h1)
        read_replace(ctl_in, name+"_alt.ctl", f1, r1)
        os.chdir(working_dir)

    print("Control files successfully created, GWidecodeml is ready for codeml performance")


    
    # step 5 : run codeml in parallel
    pool = mp.Pool(t) # number of threads
    pool.map(run_codeml, [ctl for ctl in gene_names])
    pool.close()
    os.chdir(working_dir)

    if args.lrt_stop:
        print("Codeml has finished, output files successfully created. Exiting...")
        exit()
    else:
        print("Codeml has finished, output files successfully created. Running LRTs...")


    # Step 6: LRTs, keep genes with p-value < 0.05
    alt_results = [f for f in os.listdir(working_dir) if f.endswith("_alt.txt")]
    significants = []
    for f in alt_results:
        gene_name = f.split("_")[0]
        ln_alt = lnl(f)
        ln_null = lnl(gene_name + "_null.txt")
        if model in ["BS","BM"]:
            p_val = lrt(ln_alt, ln_null)
        elif model == "SM":
            p_val = lrt(ln_alt,ln_null,2) # in site-model, LRT has 2 degrees of freedom
        if p_val < 0.05:
            significants.append(gene_name)



    # Step 7: get omega/beb values for those genes with significant values:
    out_file_name = "results_"+model+".tsv"
    out_file = open(out_file_name,"w")

    # branch model: get omega and end
    if model == "BM":
        out_file.write("Gene_name" + "\t" + "Branch_label" + "\n")
        for gene in significants:
            l_out = output_bm(gene)
            if l_out:
                out_file.write(l_out + "\n")

    # branch-site model: omega + beb positions
    elif args.mode == "BS":
        out_file.write("Gene_name" + "\t" + "AA_position" + "\t" + "AA" + "\n")
        for gene in significants:
            o =  omega_bs(gene+"_alt.txt")
            if float(o) > 1:
                l_out = output_bs(gene)
                if l_out:
                    out_file.write(l_out+"\n")


    # site model: beb positions under PS
    elif args.mode == "SM":
        out_file.write("Gene_name" + "\t" + "Nr_positions" + "\t" + "AA Positions" + "\n")
        for gene in significants:
            l_out = output_sm(gene)
            if l_out:
                out_file.write(l_out+"\n")

    out_file.close()

    print("GWideCodeML succesfully run, please check {} results file.".format(out_file_name))

if __name__ == "__main__":
    main()

