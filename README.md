# interface

#### Inspiration
This script is largely inspired from the works of Dr. Francesco Raimondi in the paper [here](https://pubmed.ncbi.nlm.nih.gov/31160049/)

#### About
For a given pair of chain IDs in a [PDB](https://www.rcsb.org/) complex, this script calculates all the interfaces between them and maps them to their equivalent positions in their corresponding proteins' FASTA sequences (UniProt Accession)

#### Requirements
*psi-blast* & *makeblastdb* from the [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/) toolkit

#### Run
```
./interface.py 3sn6 R A
```

#### Protocol
1. User enters pair of chain IDs and their PDB-ID complex
2. Fetch the PDB complex using API
3. Fetch [SIFTS](https://www.ebi.ac.uk/pdbe/api/doc/sifts.html) annotation of PDB-ID -> Chain -> UniProt accession using API
4. Calculate interfaces (considering all atoms) between the given chains
5. Map interfaces onto protein FASTA sequences using PSI-BLAST
6. Save the output

#### Arguments
1. pdb: PDB-ID of the complex
2. chain A: First chain of the PDB complex
3. chain B: Second chain of the PDB complex
4. dist: Cutoff of the interface distance (in Angstroms; default is 6.5)
5. temp: Temporary directory to store the output (user can also provide one; useful especially when the user plans to run this script multiple times on the same PDB complex)
6. out: Path & name of the output file (default: print on screen)

#### Future developments
1. Generate PyMol output
2. Select type of atoms to consider for interface (default: all types)
3. Use MMCIF format
4. Improve selection of hits in PSI-BLAST
