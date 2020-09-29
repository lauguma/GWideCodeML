## Useful scripts for your dataset

This folder contains useful scripts for creating input files for running GWideCodeML.

### cds2concat.py

This script might be useful for obtaining species phylogeny necessary for running GWideCodeML.
It takes as input a set of fasta files of orthologs. Each fasta file will contain gene sequences
corresponding to different species labeled with a conserved identifier ("SPPTAG") among all the orthologs provided.
This identifier will be used to identify orthologs having representatives of all the species provided 
in the analyses. The script will identify those orthologs and create files necessary for phylogeny reconstruction.
Sequence names must be in the following format:
'>SPPTAG_genename'

e.g.: 

`>SCERS288C_YPL091W`

`python cds2concat.py -h`

```
usage: cds2concat.py [-h] [-cds CDS_SUFFIX] [-p THREADS] [-work_dir WD]
                     [-aligner TOOL]

This script takes fasta files containing annotated genes as input: one fasta per gene
    and one species/strain per fasta.
    Fasta files extension must be indicated (-cds option) and files must be located in the working dir.
    The output will be:
    - Translated fasta files into amino-acid sequences (.faa)
    - Alignments of protein sequences with mafft (.faa.mafft)
    - Back-translated alignments into codons (.cd.mafft)
    - Alignment concatenate of containing all codon alignments (cat.cd.mafft)
    

optional arguments:
  -h, --help       show this help message and exit
  -cds CDS_SUFFIX  fasta file suffix containing cds sequences
  -p THREADS       number of threads
  -work_dir WD     working directory
  -aligner TOOL    MSA tool: prank, mafft or muscle (default)

```



### fsa2msa.py

`python fsa2msa.py -h`

```
usage: fas2msa.py [-h] [-fasta CDS_SUFFIX] [-p THREADS] [-aligner TOOL]

optional arguments:
  -h, --help         show this help message and exit
  -fasta CDS_SUFFIX  fasta file suffix containing cds sequences
  -p THREADS         number of threads
  -aligner TOOL      MSA tool: prank, mafft or muscle (default)
  ```

This script takes as input nucleotide coding sequences in fasta file format:
- Translate them into aminoacids (.faa files).  
- Align aminoacid sequences. User can choose alignment tool between prank, mafft and muscle (default). Please, note that prank is much more slower than the others (.faa.msa files).  
- Back-translate aligned aminoacid sequences into codon alignments (final output, .cd.msa files).  
