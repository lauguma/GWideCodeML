# GWideCodeML


GWideCodeML is a Python package that provides support for testing evolutionary hypothesis using codeml (from the PAML package) in a genome-wide framework.


For further information on installation and usage, please visit https://github.com/lauguma/GWideCodeML/wiki



## Requirements


Python >= 3.5

**Python libraries**

* Biopython 
* Scipy 
* ete3 

**Software**

* *codeml* from the PAML package (Phylogenetics Analysis by Maximum Likelihood). http://abacus.gene.ucl.ac.uk/software/paml.html
* fastTree 

Easy installation recommended via pip install using conda environments.

## Installation

1. Download GWideCodeML
`git clone https://github.com/lauguma/GWideCodeML.git`

2. Install
`python setup.py install`

3. Run GWideCodeML: if succesfull installation, gwidecodeml executable is created. You can check it by writing in your console:
`gwidecodeml -h`


## Citation

Macías L. G., Barrio E. and Toft. C. GWideCodeML: a Python package for testing evolutionary hypothesis at the genome-wide level (unpublished).

Our pipeline uses third-party software. If you are using GWideCodeML, please, also cite:

Yang, Z. "PAML 4: a program package for phylogenetic analysis by maximum likelihood."
Mol Biol Evol (2017) doi: 10.1093/molbev/msm088 

Huerta-Cepas, J., Serra, F and Bork, P. "ETE 3: Reconstruction,
analysis and visualization of phylogenomic data."  Mol Biol Evol (2016) doi:
10.1093/molbev/msw046


