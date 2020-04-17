#!/usr/bin/python

'''
GWideCodeML: a Python automatic pipeline for running codeml
'''

__author__ = "Laura Gutierrez Macias and Christina Toft"
__version__ = "1.0.1"
__mantainer__ = "Laura G. Macias"
__email__ = "laugmacias@gmail.com"
__status__ = "Development"

import os
import subprocess
from argparse import ArgumentParser
import multiprocessing as mp
from scipy import stats
from ete3 import EvolTree
from Bio import SeqIO
import re



# Genes must be annotated with a 4 letters tag followed by strain name + _ + gene name
# example: SUVAZP964_12G2550
def species_name(name_string):
    spcode = name_string[:4]
    return spcode

# Strain name is followed by _
def strain_name(name_string):
    #strcode = name_string.split("_")[0]
    strcode = name_string
    return strcode

# Read multi-fasta file and return sequence IDs in a list
def strain_ids(fasta):
    ids = []
    fasta_file = SeqIO.parse(fasta, "fasta")
    for s in fasta_file:
        ids.append(strain_name(s.id.split("_")[0]))
    fasta_file.close()
    return ids


# Find all strain from the same spp (4 letters TAG)
def find_strains(str_list, spp_tag):
    found = []
    for x in str_list:
        if x.startswith(spp_tag):
            found.append(x)
    return found

# Transform fasta format to format compatible with paml
def fasta2phy(msa_input, phy_out):
    input_handle = open(msa_input, "rU")
    output_handle = open(phy_out, "w")
    seqs = []
    headers = []
    alignments = SeqIO.parse(input_handle, "fasta")
    for a in alignments:
        headers.append(str(a.id))
        seqs.append(str(a.seq))
    input_handle.close()
    # write first line of the output file
    output_handle.write("  " + str(len(headers)) + "  " + str(len(seqs[0])) + "  " + "\n")
    # write sequences in paml format
    for x in range(0, len(headers)):
        output_handle.write(headers[x] + "  " + seqs[x] + "\n")
    output_handle.close()

# Extract Ln Likelihood values from output files created by codeml
# *** anhadir error handle if file is empty
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
def modelSelection(select):
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
        "omega":"1.4",
        "ncatG":"3"
    }
    # branch model null hypothesis
    bm0 = {
        "model":"0",
        "NSsites":"0",
        "fix_omega":"0",
        "omega":"1.4",
        "ncatG":"3"
    }
    # branch model alt hypothesis
    bm1 = {
        "model":"2",
        "NSsites":"0",
        "fix_omega":"0",
        "omega":"1.4",
        "ncatG":"3"
    }
    # site model M1a null hypothesis
    sm0 = {
        "model":"0",
        "NSsites":"1",
        "fix_omega":"0",
        "omega":"1.4",
        "ncatG":"3"
    }
    # site model M2a alt hypothesis
    sm1 = {
        "model":"0",
        "NSsites":"2",
        "fix_omega":"0",
        "omega":"1.4",
        "ncatG":"3"
    }

    if select == "BS":
        return bs0,bs1
    elif select == "BM":
        return bm0,bm1
    elif select == "SM":
        return sm0,sm1
    else:
        print "Error, incorrect model selected"


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
        print "Error: codeml.ctl. No such file in the working directory"


# Create and edit codeml.ctl file to run codeml according to model selected
def codeml_settings(seq_input, tree_input, out_name, model):

    ctl_find = ["seq_file.phy","tree.nwk","out_name","model = 2","NSsites = 2","fix_omega = 0","omega = .4","ncatG = 3"]
    ctl_replace = [seq_input,tree_input,out_name]

    ctl_replace.append("model = " + model["model"])
    ctl_replace.append("NSsites = " + model["NSsites"])
    ctl_replace.append("fix_omega = " + model["fix_omega"])
    ctl_replace.append("omega = " + model["omega"])
    ctl_replace.append("ncatG = " + model["ncatG"])

    return ctl_find,ctl_replace


# For codeml running, create individual folders on the working dir to run codeml on them, both null and alt
def runCodeml(ctlFile):
    os.chdir(os.path.join(working_dir,ctlFile))
    # codeml null hypothesis
    bashCommand = "codeml " + ctlFile + "_null.ctl"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()

    bashCommand = "codeml " + ctlFile + "_alt.ctl"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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

# get beb positions if there is any and their probabilities
def beb_bs(alt_file):
    f_in = open(alt_file, "r")
    pos_pv = []
    pos = []
    for line in f_in:
        if line.startswith("Bayes Empirical Bayes (BEB)"):
            f_in.next()
            line = f_in.next()
            i = True
            if line != "\n":
                while i:
                    pos_pv.append(line)  # add specific lines to matrix
                    line = f_in.next()

                    if line == "\n":
                        i = False  # stop reading when there is a \n

    if len(pos_pv) > 0:
        for x in pos_pv:
            if '*' in x:
                pos.append(int(re.search(r'\d+', x).group()))

    if len(pos) > 0:
        return pos
    else:
        return False

# get beb positions for site models, return a dictionary
def beb_sm(alt_file):
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
                    aas[text_list[0]] = text_list[2]

    return aas


# Step 1: Create database: multifasta file
# INPUT file: DB_CDS.fna


# Step 2: orthology information
# INPUT file: pillar_tab.txt


# Step 3: Read multifasta files and spp_tree
# Count species in multifasta file and prune the others
# Find common ancestor and mark branches


def main():
    parser = ArgumentParser(description="parameters for GWideCodeml performance")
    parser.add_argument("-p", dest="threads", help="number of threads", type=str,default=1)
    parser.add_argument("-model", dest="mode", help="model testing", type=str, default="BS")
    parser.add_argument("-cds", dest="suffix", help="codon alignment suffix", type=str,default=".cd.mafft")
    parser.add_argument("-branch", dest="mark", help="text file containing species of the branch of interest", type=str)
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-tree", dest="spp_tree", help="species tree in newick format", type=str)
    requiredNamed.add_argument("-work_dir", dest="wd", help="working_dir", type=str)

    args = parser.parse_args()

    # define working dir
    global working_dir
    working_dir = os.path.abspath(args.wd) # get absolute path to avoid errors

    # nr. threads define by the user
    t = int(args.threads)

    # list codon alignments in working dir
    suffix = args.suffix
    fasta_files = [f for f in os.listdir(working_dir) if f.endswith(suffix)]

    # mark branches if model selected is not site model
    # ** read file containing branches to mark
    if args.mode != "SM":
        mark_spp = "SCER"  # spp tag to mark (must be 4 letter tag)
        mark_spp2 = "SCER"

    '''
    # Parse species tree
    spptree = EvolTree(args.spp_tree)
    nodes = [] # save in a list all spp contained in spp tree
    for leaf in spptree:
        nodes.append(leaf.name)

    alignments = []
    for fasta in fasta_files:
        tree = EvolTree(args.spp_tree)  # init tree every time a fasta is open
        name = fasta.split(".")[0]
        os.mkdir(name)
        os.chdir("./" + name)




        genomes = strain_ids("../"+fasta)  # seqs in fasta

        # 1: prune tree
        spp2tree = list(set(nodes).intersection(genomes)) # spp to keep in tree
        tree.prune(spp2tree)

        # 2: mark branches
        # en -model SM no hay que marcar ramas, solo podar


        #branch = find_strains(genomes,mark_spp)
        #branch2 = find_strains(genomes,mark_spp2)
        #ancestor = tree.get_common_ancestor(branch)  # get ancestor of spp to mark branches
        #ancestor2 = tree.get_common_ancestor(branch2)
        #n1 = str(ancestor.node_id)
        #n2 = str(ancestor2.node_id)
        #tree.mark_tree([n1,n2], marks=["#1","#2"])
        #tree.mark_tree([n1], marks=["#1"])
        tree.write(outfile=name+".tree") # write tree with only topology
        #print tree.write()


        # File format converter: MSA fasta --> Phylip
        fasta2phy("../"+fasta, name+".phy")


    # step 4: create codeml.ctl files
        ctl_in = os.path.join(working_dir,"codeml.ctl")
        h0,h1 = modelSelection(args.mode)
        # create ctl file for null hypothesis testing
        f0,r0 = codeml_settings(name+".phy", name+".tree", name+"_null.txt", h0)
        read_replace(ctl_in, name+"_null.ctl", f0, r0)

        # create ctl file for alternative hypothesis testing
        f1, r1 = codeml_settings(name + ".phy", name + ".tree", name + "_alt.txt", h1)
        read_replace(ctl_in, name+"_alt.ctl", f1, r1)
        os.chdir(working_dir)
        alignments.append(name)

    print "Control files successfully created, GWidecodeml is ready for codeml performance"


    # step 5 : run codeml in parallel
    pool = mp.Pool(t) # number of threads
    pool.map(runCodeml, [ctl for ctl in alignments])
    pool.close()
    os.chdir(working_dir)
    '''
    # Step 6: get likelihood values and calculate p-value by LRT

    alt_results = [f for f in os.listdir(working_dir) if f.endswith("_alt.txt")]
    significants = []

    for f in alt_results:
        gene_name = f.split("_")[0]
        ln_alt = getLn(f)
        ln_null = getLn(gene_name + "_null.txt")
        p_val = lrt(ln_alt, ln_null)
        if p_val < 0.05:
            significants.append(gene_name)


    # Step 7: get omega/beb values for those genes with significant values:

    for gene in significants:

        # branch model: get omega and end
        if args.mode == "BM":
            o = omega_bm(gene+"_alt.txt")
            if float(o[2]) > 1:
                print gene


        # branch-site model: omega + beb positions
        elif args.mode == "BS":
            o =  omega_bs(gene+"_alt.txt")
            sites = beb_bs(gene+"_alt.txt")
            if float(o) > 1:
                print o,sites

        #site: omega + beb positions
        elif args.mode == "SM":
            print gene
            positions = beb_sm(gene+"_alt.txt")
            print positions



if __name__ == "__main__":
    main()

