#!/usr/bin/python

import argparse
from Bio import SeqIO
from Bio import Seq
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
import multiprocessing as mp
import subprocess
import os



# USAGE python ffn2align -h

def read_spptag(fasta_in):
    """Read a fasta file a returns species names contained in headers
    Species name format: >SPPTAG_genename e.g.: >SCERS288C_YPL091W"""
    tags = []
    fasta_handle = SeqIO.parse(fasta_in, "fasta")
    for s in fasta_handle:
        # Change here if fasta header is in a different format #
        tags.append(s.id.split("_")[0])
    fasta_handle.close()
    return tags

def remove_duplicates(list_dup):
    """Remove duplicated elements in list"""
    list_out = list(dict.fromkeys(list_dup))
    return list_out

def list_files(work_dir, suffix):
    """Return a list of files with a certain pattern in a directory"""
    file_list = [f for f in os.listdir(work_dir) if f.endswith(suffix)]
    return file_list

def seqs_in_fasta(seqs_list,fasta_in):
    """Returns true if all sequence names from a list are in fasta"""
    fasta_names = read_spptag(fasta_in)
    if set(fasta_names) == set(seqs_list):
        return True
    else:
        return False

def translate(nts_in):
    """This function opens a fasta file containing cds, read,
    translate and create another fasta file with out_file as name"""
    out_file = nts_in.replace(".ffn", ".faa")
    fasta = SeqIO.parse(nts_in, "fasta")
    prot = open(out_file, "w")
    for s in fasta:
        aas = s.seq.translate(stop_symbol="")
        sequence_object = Seq(str(aas), IUPAC.ExtendedIUPACProtein)
        record = SeqRecord(sequence_object, id=s.id, description="")
        SeqIO.write(record, prot, "fasta")
    fasta.close()
    prot.close()


# run mafft: output aligned file written to the standard output
def run_mafft(aa_file):
    bash_command = "mafft --quiet " + aa_file
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    process_output, error = process.communicate()
    file = open(aa_file+".msa", "w")
    file.write(process_output.decode('utf-8'))
    file.close()

# run muscle
def run_muscle(aa_file):
    bash_command = "muscle -in " + aa_file + " -out " + aa_file + ".msa"
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    output = process.stdout.read()

# run prank: (slower than mafft and muscle)
def run_prank(aa_file):
    gname = aa_file.replace(".faa","")
    bash_command = "prank -d=" + aa_file + " -o=" + gname
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    output = process.stdout.read()

    os.rename(gname+".best.fas", gname+".faa.msa")


def back_translate(ali_file):
    """Read aligned amino acid sequences and back translate
    into codon alignments by reading DNA sequences"""
    gene_name = ali_file.replace(".faa.msa", "")
    ali = SeqIO.parse(ali_file, "fasta")
    out_file = open(gene_name + ".cd.msa", "w")
    for s in ali:
        aa_seq = str(s.seq)
        dna = SeqIO.parse(gene_name + gene_suffix, "fasta")
        for cds in dna:
            if cds.id == s.id:
                # print cds.id
                nt_seq = str(cds.seq)
                c1 = 0  # codon start
                c2 = 3  # codon end
                chain = ""
                for aa in aa_seq:
                    if aa != "-":
                        codon = nt_seq[c1:c2]
                        c1 += 3
                        c2 += 3
                        chain += codon
                    else:
                        chain += aa * 3

                sequence_object = Seq(chain, IUPAC.unambiguous_dna)
                record = SeqRecord(sequence_object, id=cds.id, description="")
                SeqIO.write(record, out_file, "fasta")
        dna.close()


def cat_msa(suffix_msa):
    """Concatenate back-translated aligned cds sequences"""
    path = os.getcwd()  # working directory
    multi_files = [f for f in os.listdir(path) if f.endswith(suffix_msa)]

    # create dict with all sequences to be concatenated
    spp = dict()
    spp_names = read_spptag(multi_files[0])
    for x in spp_names:
        spp[x] = ""
    for f in multi_files:
        fasta = SeqIO.parse(f, "fasta")
        for s in fasta:
            # change here if fasta header format is different #
            tag = s.id.split("_")[0]
            spp[tag] += str(s.seq)
        fasta.close()
    cat_file = open("cat.cd.msa", "w")
    for k, v in spp.items():
        sequence_object = Seq(str(v), IUPAC.unambiguous_dna)
        record = SeqRecord(sequence_object, id=k, description="")
        SeqIO.write(record, cat_file, "fasta")
    cat_file.close()


def main():
    """This script takes fasta files containing annotated genes as input: one fasta per gene
    and one species/strain per fasta, all the fasta files must contain the same species/strains
    named under the same way.
    Fasta files extension must be indicated (-cds option) and files must be located in the working dir.
    The output will be:
    - Translated fasta files into amino-acid sequences (.faa)
    - Alignments of protein sequences with mafft (.faa.mafft)
    - Back-translated alignments into codons (.cd.mafft)
    - Alignment concatenate of containing all codon alignments (cat.cd.mafft)
    """

    parser = argparse.ArgumentParser(description=main.__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-cds", dest="cds_suffix", help="fasta file suffix containing cds sequences",
                        type=str, default=".ffn")
    parser.add_argument("-p", dest="threads", help="number of threads", type=str, default=1)
    parser.add_argument("-work_dir", dest="wd", help="working directory", type=str,
                        default=os.path.abspath(os.getcwd()))
    parser.add_argument("-aligner", dest="tool", help="MSA tool: prank, mafft or muscle (default)",
                        type=str, default="muscle")
    args = parser.parse_args()

    t = int(args.threads)
    path = os.path.abspath(args.wd)   # working directory
    global gene_suffix
    gene_suffix = args.cds_suffix

    # Step 0: identify fasta file containing orthologs for all species
    nt_files = list_files(path, gene_suffix)
    species = []
    for x in nt_files:
        species.extend(read_spptag(x))
    all_spp = remove_duplicates(species)

    genes_to_cat = []
    for x in nt_files:
        if seqs_in_fasta(all_spp,x):
            genes_to_cat.append(x)
    print("Number of genes having all species names: {}".format(str(len(genes_to_cat))))
    # Step 1: translate nt (suffix) to aa (.faa)
    pool = mp.Pool(t)  # number of threads
    pool.map(translate, [ffn for ffn in genes_to_cat])
    pool.close()

    # Step 2: align protein sequences
    aa_files = list_files(path, ".faa")
    pool = mp.Pool(t)  # number of threads
    if args.tool == "muscle":
        pool.map(run_muscle, [faa for faa in aa_files])
        pool.close()
    elif args.tool == "mafft":
        pool.map(run_mafft, [faa for faa in aa_files])
        pool.close()
    elif args.tool == "prank":
        pool.map(run_prank, [faa for faa in aa_files])
        pool.close()
    else:
        print("Wrong aligner selected. Please, choose between mafft, prank or muscle (default)")
        exit()

    # Step 3 : back-translate into codons

    ali_files = list_files(path, ".faa.msa")
    pool = mp.Pool(t)  # number of threads
    pool.map(back_translate, [msa for msa in ali_files])
    pool.close()

    # Step 4: concatenate codon-alignments
    cat_msa(".cd.msa")
    # OUTPUT: cat.cd.msa --> concatenate alignment in fasta format ready to use for phylogeny (e.g. RAxML)

if __name__ == "__main__":
    main()
