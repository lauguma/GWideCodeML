#!/usr/bin/python

from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import multiprocessing as mp
import subprocess
import os


# USAGE python fas2msa.py -h

'''
This script takes fasta files containing annotated genes as input: one fasta per gene 
and one species/strain per fasta, all the fasta files must contain the same species/strains 
named in the same way. 
Fasta files extension must be indicated (-cds option) and files must be located in the working dir.
The output will be:
- Translated fasta files into amino-acid sequences (.faa)
- Alignments of protein sequences with mafft (.faa.mafft)
- Backtranslated alignments into codons (.cd.mafft)
'''


# function translate: opens a fasta file containing cds, read,
# translate and create another fasta file with out_file as name
def translate(nts_in):
    suffix_fasta = nts_in.split(".")[-1]
    out_file = nts_in.replace(suffix_fasta, "faa")
    fasta = SeqIO.parse(nts_in, "fasta")
    prot = open(out_file, "w")

    for s in fasta:
        aas = s.seq.translate(stop_symbol="")
        sequence_object = Seq(str(aas))
        record = SeqRecord(sequence_object, id=s.id, description="")
        SeqIO.write(record, prot, "fasta")

    prot.close()


# run mafft: output aligned file written to the standard output
def run_mafft(aa_file):
    bashCommand = "mafft-linsi " + aa_file
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output = process.stdout.read()

    with open(aa_file + ".msa", "w") as f_out:
        f_out.write(output.decode('utf-8'))

# run muscle
def run_muscle(aa_file):
    bashCommand = "muscle -in " + aa_file + " -out " + aa_file + ".msa"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output = process.stdout.read()

# run prank: (slower than mafft and muscle)
def run_prank(aa_file):
    gname = aa_file.replace(".faa","")
    bashCommand = "prank -d=" + aa_file + " -o=" + gname
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output = process.stdout.read()

    os.rename(gname+".best.fas", gname+".faa.msa")


def backTranslate(gene_file):
    gene_name = gene_file.split(".")[0]
    msa = SeqIO.parse(gene_name+".faa.msa", "fasta")
    out_file = open(gene_name + ".cd.msa", "w")

    for s in msa:
        aa_seq = str(s.seq)
        dna = SeqIO.parse(gene_file, "fasta")
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

                sequence_object = Seq(chain)
                record = SeqRecord(sequence_object, id=cds.id, description="")
                SeqIO.write(record, out_file, "fasta")


def main():
    parser = ArgumentParser()
    parser.add_argument("-fasta", dest="cds_suffix", help="fasta file suffix containing cds sequences", type=str)
    parser.add_argument("-p", dest="threads", help="number of threads", type=str,default=1)
    parser.add_argument("-aligner", dest="tool", help="MSA tool: prank, mafft or muscle (default)", \
                        type=str, default="muscle")

    args = parser.parse_args()

    t = int(args.threads)
    path = os.path.abspath(os.getcwd())  # working directory

    # gene fasta file suffix
    genes_suffix = args.cds_suffix


    # Step 1: translate nt (.ffn) to aa (.faa)
    nt_files = [f for f in os.listdir(path) if f.endswith(genes_suffix)]

    pool = mp.Pool(t)  # number of threads
    pool.map(translate, [ffn for ffn in nt_files])
    pool.close()

    # Step 2: align protein sequences
    aa_files = [f for f in os.listdir(path) if f.endswith(".faa")]

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




    # Step 3 : backtranslate into codons

    pool = mp.Pool(t)  # number of threads
    pool.map(backTranslate, [x for x in nt_files])
    pool.close()

if __name__ == "__main__":
    main()
