#!/usr/bin/python

import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import shutil



def genes2list(fasta_file):
    genes = []
    fasta = SeqIO.parse(fasta_file,"fasta")
    for s in fasta:
        genes.append(s.id)
    fasta.close()
    return genes

def parseFasta(fasta_in,fasta_out):

    fasta = SeqIO.parse(fasta_in,"fasta")
    f_out = open(fasta_out,"a")
    for s in fasta:
        sequence_object = Seq(str(s.seq), IUPAC.unambiguous_dna)
        record = SeqRecord(sequence_object, id=s.id, description="")
        SeqIO.write(record, f_out, "fasta")

    fasta.close()
    f_out.close()


def fasta2dict(fasta_in):
    fasta_dict = dict()
    fasta = SeqIO.parse(fasta_in, "fasta")
    for s in fasta:
        gene = s.id
        if ";" in s.id:
            gene = s.id.split(";")[0]

        fasta_dict[gene] = str(s.seq)

    fasta.close()
    return fasta_dict


def dict2ffn(fdict, gene_list,out_name):
    fasta_out = open(out_name+".ffn","w")

    for gene in gene_list:
        if gene in fdict.keys():
            bases = fdict[gene]
            sequence_object = Seq(bases, IUPAC.unambiguous_dna)
            record = SeqRecord(sequence_object, id=gene, description="")
            SeqIO.write(record, fasta_out, "fasta")

    fasta_out.close()


def listReplace(nlist,element,replacement):
    for n, i in enumerate(nlist):
        if i == element:
            nlist[n] = replacement
    return nlist


def main():

    # INPUTs
    fasta_db = "DB_CDS.fna"
    fasta_db_out = "DB_new.fna"
    pillar = pd.read_csv("PILLAR.tab", delimiter="\t")

    ## 1
    # Create fasta database: ygob database + fasta.cds.fna
    path = os.getcwd() # working directory
    suffix = ".cds.fna" # suffix of the multifasta files containing annotated genomes
    fasta_files = [f for f in os.listdir(path) if f.endswith(suffix)]
    db_out = open(fasta_db_out,"w") # iniciar db para que empiece en blanco
    db_out.close()
    parseFasta(fasta_db,fasta_db_out)
    for f in fasta_files:
        parseFasta(f,fasta_db_out)


    ## 2
    # Read pillar tab and add gene names of genomes  to analyze
    preWGD = ["ZROU", "TDEL", "KLAC", "EGOS", "ECYM", "LKLU", "LTHE", "LWAL"]
    postWGD = ["SCERREF", "SMIKIFO1815", "SUVACBS7001", "CGLA", "KAFR", \
                 "KNAG", "NCAS", "NDAI", "TBLA", "TPHA", "VPOL"]
    spp1 = ["SCERT73","SCERREF","SCERY9","SCERYPS128"] # definir spp de interes
    pillar['SCERREF.1'] = pillar['SCERREF.1'].str.replace('SCERREF_', '')
    pillar['SCERREF.2'] = pillar['SCERREF.2'].str.replace('SCERREF_', '')

    pillar.fillna(0, inplace=True)


    # Separate pre and post-WGD species
    pillar1 = pillar.iloc[:,0:20]
    pillar2 = pillar.iloc[:,11:]

    pillar1.columns = pillar1.columns.str.replace('\.1', '')
    pillar2.columns = pillar2.columns.str.replace('\.2', '')

    p1 = pillar1[(pillar1[postWGD] != 0).any(1)]
    p2 = pillar2[(pillar2[postWGD] != 0).any(1)]

    mytab = pd.concat([p1,p2], axis=0, join="inner", ignore_index=True)

    # remove all-zeros rows
    mytab = mytab[(mytab != 0).any(1)]


    # gene orthologs will be added by systematic gene name ortholog to Scerevisiae S288c
    s288c = list(filter(lambda num: num != 0, mytab["SCERREF"].tolist()))

    d = dict() # dict with key strain name and list of genes as values
    for f in fasta_files:
        strain = f.replace(suffix,"")
        postWGD.append(strain)
        d["%s" %strain] = genes2list(f)

    # iterate over de dict and start filling the table
    for k,v in d.items():
        mytab[k] = 0 # create new zeros column with strain name as column name

        # iterate over gene list of every genome
        for gene in v:
            #print(gene)
            ortho = gene.split("_")[1] # ortho value to find gene in table
            if ortho in s288c:
                row = int(mytab[mytab['SCERREF'] == ortho].index[0]) # get row of ortho gene
                mytab.loc[row, k] = gene

    # save table with the new orthologs added
    mytab.to_csv("ortho_table.txt", index=False, sep="\t", encoding='utf-8')  # para comprobar

    # Read table: uncomment previous section to start from this part
    subdata = pd.read_csv("ortho_table.txt",sep = "\t",encoding='utf-8')
    cds_db = fasta2dict("DB_new.fna")

    ## 4
    # Create fasta file gene by gene
    for i in range(len(subdata.index)):
        gname = subdata.loc[i, "SCERREF"]
        r = subdata.iloc[i].tolist()
        n = listReplace(r,gname,"SCERREF_"+gname)
        dict2ffn(cds_db,n,gname)

    # Divide gene fasta file in 2 sub-folders:
    # spp_tree/ orthologs shared among all species, those genes can be used for generate the species tree
    os.mkdir("spp_tree")
    os.mkdir("other_genes")
    allspp = subdata[(subdata != "0").all(1)]
    genes_phylo = allspp["SCERREF"].tolist()
    for ffn in genes_phylo:
        ffn_file = ffn+".ffn"
        shutil.move(ffn_file,"./spp_tree/"+ffn_file)

    # other_genes/ differential orthologs among different species
    for ffn in [f for f in os.listdir(os.getcwd()) if f.endswith(".ffn")]:
        shutil.move(ffn, "./other_genes/" + ffn)



if __name__ == "__main__":
    main()
