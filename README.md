# interface

#### Inspiration
This script is largely inspired from the works of Dr. Francesco Raimondi in the *[Cell](https://pubmed.ncbi.nlm.nih.gov/31160049/)* paper

#### About
For a given pair of chain IDs in a [PDB](https://www.rcsb.org/) complex, this script calculates all the interfaces between them and maps them to their equivalent positions in their corresponding proteins' FASTA (UniProt Accession)

#### Requirements
*psi-blast* & *makeblastdb* from the [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/) toolkit

#### Run
```
./interface.py 3sn6 R A
```
#### Future developments
1. Generate PyMol output
2. Select the type of atoms to consider for interface (default: all)
