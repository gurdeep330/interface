#!/usr/bin/python3.8
'''intermapper: determine interface residues of any two chains in a protein complex'''

'''
This script is largely inspired from the works of Francesco Raimondi in
determining GPCR/G-protein interfaces.
intermapper finds all possible interface residues of any two chains in a PDB
complex. Then it maps the interface residues to their FASTA sequences in UniProt.
'''

## Import libraries
import os, sys, gzip, argparse, requests, tempfile
import xml.etree.ElementTree as ET
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import *
import numpy as np

__author__ = "Gurdeep Singh"
__credits__ = ["Gurdeep Singh", "Francesco Raimondi"]
__license__ = "GPL"
__version__ = "1.0"
__email__ = "gurdeep330@gmail.com"

## Parse the input
path = os.getcwd()
parser = argparse.ArgumentParser(description='Finds interface residues of any two '+
                                'chains in a protein complex',
                                epilog='Contact: gurdeep330[at]gmail[dot]com')
parser.add_argument('pdb', help='PDB-ID (pdb/cif format)')
parser.add_argument('chainA', help='ChainA-ID')
parser.add_argument('chainB', help='ChainB-ID')
parser.add_argument('--dist', help='interface distance cutoff (in Angstroms; '+
                                'default: 6.5)')
parser.add_argument('--temp', help='Path of directory where files created'+
                    'during the run should be stored (default is temporary directory)')
parser.add_argument('--out', help='Path of file where the output should'+
                    'be stored (default is printed on the screen)')
args = parser.parse_args()

pdb = args.pdb
given_chainA = args.chainA
given_chainB = args.chainB
DISTCUTOFF = 6.5
if args.dist != None:
    DISTCUTOFF = float(args.dist)
if args.temp != None:
    tempdir = args.temp + '/'
else:
    tempdir = tempfile.TemporaryDirectory().name
    #tempdir = tempfile.mkdtemp() + '/'

out = None
if args.out != None:
    out = args.out

## Define the class
class interface:
    def __init__(self, residueA, residueB):
        self.residueA = residueA
        self.residueB = [residueB]
        self.fasta_positionA = {}
        self.fasta_positionB = {}

    def add_residueB(self, residueB):
        self.residueB.append(residueB)

## Download the protein complex
print (f'Fetching {pdb}')
try:
    os.system('wget https://files.rcsb.org/view/'+pdb+'.cif -nv -O '+tempdir+pdb+'.cif')
except:
    os.system('wget https://files.rcsb.org/view/'+pdb+'.pdb -nv -O '+tempdir+pdb+'.pdb')

format = 'cif'
if os.path.getsize(tempdir+pdb+'.cif') == 0:
    format = 'pdb'
    if os.path.getsize(tempdir+pdb+'.pdb') == 0:
        print (f'{pdb} not found')
        sys.exit()

## 3 letters to 1 letter dictionary of Amino Acids
AA ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q',
'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',
'GLY':'G', 'PRO':'P', 'CYS':'C'}

## Funcition to map chain positions to FASTA
def map_chain_to_fasta(chain, fasta, dic):
    open(tempdir+pdb+'_'+chain.id+'.fasta', 'w').write(fasta)
    os.system('makeblastdb -dbtype prot -out '+tempdir+pdb+'_'+chain.id+' -in '+tempdir+pdb+'_'+chain.id+'.fasta')
    pdb_to_unp = {}
    for acc in dic[chain.id]:
        if os.path.isfile(tempdir+acc+'.fasta') == False:
            #os.system('wget https://www.uniprot.org/uniprot/'+acc+'.fasta')
            #os.system('mv '+acc+'.fasta '+tempdir)
            requestURL = "https://www.ebi.ac.uk/proteins/api/proteins/"+acc
            r = requests.get(requestURL, headers={ "Accept" : "application/json"})
            if not r.ok:
              r.raise_for_status()
              sys.exit()

            #print (r.json()['sequence']['sequence'])
            open(tempdir+acc+'.fasta', 'w').write(r.json()['sequence']['sequence'])
            #sys.exit()
        os.system('psiblast -out '+tempdir+acc+'_blastp.txt -query '+tempdir+acc+'.fasta -db '+tempdir+pdb+'_'+chain.id)
        pdb_to_unp[acc] = {}
        for line in open(tempdir+acc+'_blastp.txt', 'r'):
            if len(line.split()) > 0:
                if line.split()[0] == 'Query':
                    sq = int(line.split()[1])
                    query = line.split()[2]
                    eq = int(line.replace('\n', '').split()[-1])
                elif line.split()[0] == 'Sbjct':
                    ss = int(line.split()[1])
                    sbjct = line.split()[2]
                    es = int(line.replace('\n', '').split()[-1])

                    for q, s in zip(query, sbjct):
                        if s not in ['.', '-'] and q not in ['.', '-']:
                            if ss not in pdb_to_unp[acc]:
                                pdb_to_unp[acc][ss] = sq
                            ss += 1
                            sq += 1
                        elif q in ['.', '-']:
                            ss += 1
                        elif s in ['.', '-']:
                            sq += 1
        #print (unp_to_pfam['P29274'][292])
    return pdb_to_unp

## Call the main funcition
def main():
    ## Store PDB ->  Chain -> UniProt Acc in a dictionary extracted from SIFTS
    print ('Fetching the SIFTS annotations of PDB -> Chain -> UniProt accession')
    dic = {}
    r = requests.get(url="https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/"+pdb)
    #print(r.json()['3sn6']['UniProt'])
    for acc in r.json()[pdb]['UniProt']:
        for char in r.json()[pdb]['UniProt'][acc]['mappings']:
            #print (acc, char['chain_id'])
            chain = char['chain_id']
            if chain not in dic:
                dic[chain] = []
            if acc not in dic[chain]:
                dic[chain].append(acc)

    print ('Calculating interfaces with cutoff', str(DISTCUTOFF), '...')

    ## To check if both chains are present in the complex
    if format == 'cif':
        try:
            parser = MMCIFParser()
            structure = parser.get_structure(pdb, tempdir+pdb+'.cif')
        except:
            os.system('wget https://files.rcsb.org/view/'+pdb+'.pdb -nv -O '+tempdir+pdb+'.pdb')
            parser = PDBParser()
            structure = parser.get_structure(pdb, tempdir+pdb+'.pdb')
    else:
        parser = PDBParser()
        structure = parser.get_structure(pdb, tempdir+pdb+'.pdb')

    for model in structure:
        if given_chainA in model and given_chainB in model:
            chainA = model[given_chainA]
            chainB = model[given_chainB]
        elif given_chainA not in model:
            print (f'Chain {given_chainA} absent in {pdb}')
            sys.exit()
        elif given_chainB not in model:
            print (f'Chain {given_chainB} absent in {pdb}')
            sys.exit()
        else:
            print (f'Chains {given_chainA} and {given_chainB} absent in {pdb}')
            sys.exit()

        fastaA = '>'+pdb+'_'+chainA.id+'\n'
        fastaB = '>'+pdb+'_'+chainB.id+'\n'
        map_fasta = {chainA.id: {}, chainB.id: {}}

        dic_interface = {}
        count = 0
        numA = 0
        numB = 0
        for count, residueA in enumerate(chainA):
            flag = 0
            for residueB in chainB:
                if residueA.get_resname() in AA and residueB.get_resname() in AA:
                    for atomA in residueA:
                        for atomB in residueB:
                            diff_VdW  = atomA.coord - atomB.coord
                            if np.sqrt(np.sum(diff_VdW * diff_VdW)) <= DISTCUTOFF:
                                flag = 1
                                if residueA.id[1] not in dic_interface:
                                    dic_interface[residueA.id[1]] = interface(residueA.id[1], residueB.id[1])
                                elif residueB.id[1] not in dic_interface[residueA.id[1]].residueB:
                                    dic_interface[residueA.id[1]].add_residueB(residueB.id[1])
                                break
                        if flag == 1:
                            flag = 0
                            break

            if residueA.id[0] == ' ':
                #print (pdb, chain.id, residue.get_resname(), residue.id)
                fastaA += protein_letters_3to1[residueA.get_resname()]
                numA += 1
                map_fasta[chainA.id][residueA.id[1]] = numA

        for residueB in chainB:
            if residueB.id[0] == ' ':
                #print (pdb, chain.id, residue.get_resname(), residue.id)
                fastaB += protein_letters_3to1[residueB.get_resname()]
                numB += 1
                map_fasta[chainB.id][residueB.id[1]] = numB
        #break

        #print (dic_interface)
        #sys.exit()
        pdb_to_unpA = map_chain_to_fasta(chainA, fastaA, dic)
        pdb_to_unpB = map_chain_to_fasta(chainB, fastaB, dic)
        #print ('PDB/Interface_A\tPDB_FASTA_A\tUniProt_ACC_A/Pos\tPDB/Interface_B\tPDB_FASTA_B\tUniProt_ACC_B/Pos')
        l = '# interfaces (v1.0)'
        l += '# contact: gurdeep330[at]gmail[dot]com'
        l = '# X: Chain-ID\n'
        l += '# Interface_X: Position in the PDB file\n'
        l += '# PDB_FASTA_X: Position in the PDB-FASTA file (created using C-alpha atoms)\n'
        l += '# UniProt_ACC_X: Position in the protein FASTA sequence\n'
        l += '# Interface_A\tPDB_FASTA_A\tUniProt_ACC_A\tInterface_B\tPDB_FASTA_B\tUniProt_ACC_B\n'
        for accA in dic[chainA.id]:
            for residueA in dic_interface:
                pdb_positionA = map_fasta[chainA.id][residueA]
                #pdb_positionA = dic_interface[chainA.id][residueA].fasta_positionA[residueA]
                if pdb_positionA in pdb_to_unpA[accA] and pdb_positionA != '':
                    unp_positionA = pdb_to_unpA[accA][pdb_positionA]
                    #print (dic_interface[chainA.id][residueA].residueB)
                    for accB in dic[chainB.id]:
                        for residueB in dic_interface[residueA].residueB:
                            pdb_positionB = map_fasta[chainB.id][residueB]
                            if pdb_positionB in pdb_to_unpB[accB] and pdb_positionB != '':
                                unp_positionB = pdb_to_unpB[accB][pdb_positionB]
                                #print (pdb+'/'+str(residueA)+'\t'+str(pdb_positionA)+'\t'+accA+'/'+str(unp_positionA)+'\t'+pdb+'/'+str(residueB)+'\t'+str(pdb_positionB)+'\t'+accB+'/'+str(unp_positionB))
                                l += pdb+'/'+chainA.id+'/'+str(residueA)+'\t'+str(pdb_positionA)+'\t'+accA+'/'+str(unp_positionA) + '\t'
                                l += pdb+'/'+chainB.id+'/'+str(residueB)+'\t'+str(pdb_positionB)+'\t'+accB+'/'+str(unp_positionB) + '\n'
                            else:
                                l += pdb+'/'+chainA.id+'/'+str(residueA)+'\t'+str(pdb_positionA)+'\t'+accA+'/'+str(unp_positionA) + '\t'
                                l += pdb+'/'+chainB.id+'/'+str(residueB)+'\t'+str('-')+'\t'+accB+'/'+str('-') + '\n'

                else:
                    for accB in dic[chainB.id]:
                        for residueB in dic_interface[residueA].residueB:
                            pdb_positionB = map_fasta[chainB.id][residueB]
                            if pdb_positionB in pdb_to_unpB[accB] and pdb_positionB != '':
                                unp_positionB = pdb_to_unpB[accB][pdb_positionB]
                                #print (pdb+'/'+str(residueA)+'\t'+str(pdb_positionA)+'\t'+accA+'/'+str(unp_positionA)+'\t'+pdb+'/'+str(residueB)+'\t'+str(pdb_positionB)+'\t'+accB+'/'+str(unp_positionB))
                                l += pdb+'/'+chainA.id+'/'+str(residueA)+'\t'+str('-')+'\t'+accA+'/'+str('-') + '\t'
                                l += pdb+'/'+chainB.id+'/'+str(residueB)+'\t'+str(pdb_positionB)+'\t'+accB+'/'+str(unp_positionB) + '\n'
                            else:
                                l += pdb+'/'+chainA.id+'/'+str(residueA)+'\t'+str('-')+'\t'+accA+'/'+str('-') + '\t'
                                l += pdb+'/'+chainB.id+'/'+str(residueB)+'\t'+str('-')+'\t'+accB+'/'+str('-') + '\n'

    print ('Calculation complete')
    if out != None:
        open(out, 'w').write(l)
        print ('Output written at', out)
    else:
        print (l)

main()
