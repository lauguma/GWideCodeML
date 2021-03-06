# GWideCodeML


GWideCodeML is a Python package that provides support for testing evolutionary hypothesis using codeml (from the PAML package) in a genome-wide framework.


For further information on installation and usage, please visit https://github.com/lauguma/GWideCodeML/wiki

## Installation
### Option 1:

1. Download GWideCodeML \
`git clone https://github.com/lauguma/GWideCodeML.git` \
`cd GWideCodeML`

2. Install \
`python setup.py install`

3. Run GWideCodeML: if succesfull installation, gwidecodeml executable is created. You can check it by writing in your console: \
`gwidecodeml -h` 

Note: if you choose this option, all requirements must be satisfaied before running GWideCodeML (e.g. codeml must be installed and available from the working directory), see Requirements section.

### Option 2: conda environment
(easier, recommended option)

0. Install and initialize miniconda \
(skip in case you already have a conda env and know how it works) \
https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html \
close and open a new console before continue 

1. Download GWideCodeML \
`git clone https://github.com/lauguma/GWideCodeML.git` \
`cd GWideCodeML`

2. Create a conda environment from yml file \
`conda env create -f gwidecodeml_conda.yml` 

3. Activate conda environment \
`conda activate gwcodeml`

4. Install and run GWideCodeML \
`python setup.py install` \
`gwidecodeml -h` 

## Requirements


Python >= 3.5

**Python libraries**

* Biopython 
* Scipy 
* ete3 

**Software**

* *codeml* from the PAML package (Phylogenetics Analysis by Maximum Likelihood). http://abacus.gene.ucl.ac.uk/software/paml.html \
* fastTree: http://www.microbesonline.org/fasttree/


## Citation

Macías L. G., Barrio E. and Toft. C. "GWideCodeML: a Python package for testing evolutionary hypothesis at the genome-wide level" G3: Genes, Genomes, Genetics (2020) doi:10.1534/g3.120.401874. 

Our pipeline uses third-party software:

Yang, Z. "PAML 4: a program package for phylogenetic analysis by maximum likelihood."
Mol Biol Evol (2017) doi: 10.1093/molbev/msm088 

Huerta-Cepas, J., Serra, F and Bork, P. "ETE 3: Reconstruction,
analysis and visualization of phylogenomic data."  Mol Biol Evol (2016) doi:
10.1093/molbev/msw046


