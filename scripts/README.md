### fsa2msa.py

python fsa2msa.py -h

``
usage: fas2msa.py [-h] [-fasta CDS_SUFFIX] [-p THREADS] [-aligner TOOL]

optional arguments:
  -h, --help         show this help message and exit
  -fasta CDS_SUFFIX  fasta file suffix containing cds sequences
  -p THREADS         number of threads
  -aligner TOOL      MSA tool: prank, mafft or muscle (default)
  ``

This script takes as input nucleotide coding sequences in fasta file format:
- Translate them into aminoacids (.faa files).
- Align aminoacid sequences. User can choose alignment tool between prank, mafft and muscle (default). Please, note that prank is much more slower than the others (.faa.msa files).
- Back-translate alined aminoacid sequences into codon alignments (final output, .cd.msa files).
