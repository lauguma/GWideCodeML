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
    fasta_db_out = "DB_ALL.fna"
    pillar = pd.read_csv("PILLAR.tab", delimiter="\t")

    ## 1
    # Create fasta database: ygob database + fasta.cds.fna
    path = os.getcwd() # working directory

    suffix = ".cds.fna"
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



    # Separar tablas pre-post duplicacion y despues pegar
    pillar1 = pillar.iloc[:,0:20]
    pillar2 = pillar.iloc[:,11:]

    pillar1.columns = pillar1.columns.str.replace('\.1', '')
    pillar2.columns = pillar2.columns.str.replace('\.2', '')

    p1 = pillar1[(pillar1[postWGD] != 0).any(1)]
    p2 = pillar2[(pillar2[postWGD] != 0).any(1)]

    mytab = pd.concat([p1,p2], axis=0, join="inner", ignore_index=True)

    # remove all-zeros rows
    mytab = mytab[(mytab != 0).any(1)]


    # ortologos seran anadidos por ortologo con s288c
    # anotacion sistematica
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
            #print gene
            ortho = gene.split("_")[1] # ortho value to find gene in table
            if ortho in s288c: # quitar esto y poner solo genes con ortologos con s288c, no tiene mucho sentido
                row = int(mytab[mytab['SCERREF'] == ortho].index[0]) # get row of ortho gene
                mytab.loc[row, k] = gene

    '''
    # mytab is complete with new genomes, now filter

    # iterate over rows:
    rows = [] # indexes of rows to keep in data frame

    ## 3 Este paso se hara dentro del pipeline GWideCodeml
    # Subset data by nr of outgroups and nr of genes of spp of interest
    for i in range(0, len(mytab.index)):
        # Step 1: filter by nr of outgroups
        x = mytab[preWGD].iloc[i].tolist()
        y = mytab[spp1].iloc[i].tolist()
        if y.count(0) <= 1:
            if x.count(0) <= 2:
                rows.append(i)
    
    subdata = mytab.iloc[rows,:]
    subdata.to_csv("prueba.txt",index=False,sep = "\t",encoding='utf-8') # para comprobar
    '''
    # guardar tabla completa con orthologs anadidos
    mytab.to_csv("ortho_table.txt", index=False, sep="\t", encoding='utf-8')  # para comprobar

    # Cargar base de datos en memoria
    subdata = pd.read_csv("ortho_table.txt",sep = "\t",encoding='utf-8')

    cds_db = fasta2dict("DB_ALL.fna")




    ## 4
    # Create gene multi-fasta files
    for i in range(len(subdata.index)):
        gname = subdata.loc[i, "SCERREF"]
        r = subdata.iloc[i].tolist()
        n = listReplace(r,gname,"SCERREF_"+gname)
        dict2ffn(cds_db,n,gname)



    os.mkdir("spp_tree")
    os.mkdir("other_genes")
    allspp = subdata[(subdata != "0").all(1)]
    genes_phylo = allspp["SCERREF"].tolist()
    for ffn in genes_phylo:
        ffn_file = ffn+".ffn"
        shutil.move(ffn_file,"./spp_tree/"+ffn_file)

    # mover el resto a otra carpeta
    for ffn in [f for f in os.listdir(os.getcwd()) if f.endswith(".ffn")]:
        shutil.move(ffn, "./other_genes/" + ffn)



if __name__ == "__main__":
    main()
