import argparse
import sys
import os
import numpy as np
import pandas as pd
import glob
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
import pdb

parser = argparse.ArgumentParser(description = '''Read, clean and extract sequence from a protein from a pdb file.''')

parser.add_argument('--structure_file', nargs=1, type= str, default=sys.stdin, help = 'Path to protein structure (PDB).')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')

##############FUNCTIONS##############
def read_pdb(pdbname):
    '''Read PDB
    '''

    three_to_one = {'ARG':'R', 'HIS':'H', 'LYS':'K', 'ASP':'D', 'GLU':'E',
                    'SER':'S', 'THR':'T', 'ASN':'N', 'GLN':'Q', 'CYS':'C',
                    'GLY':'G', 'PRO':'P', 'ALA':'A', 'ILE':'I', 'LEU':'L',
                    'MET':'M', 'PHE':'F', 'TRP':'W', 'TYR':'Y', 'VAL':'V',
                    'SEC':'U', 'PYL':'O', 'GLX':'X', 'UNK': 'X'}

    f=open(pdbname,'rt')

    if '.pdb' in pdbname:
        parser = PDBParser()
        struc = parser.get_structure('', f)
    else:
        parser = MMCIFParser()
        struc = parser.get_structure('',f)

    #Save
    model_seqs = {}


    for model in struc:
        for chain in model:
            #Save

            model_seqs[chain.id]=''


            #Go through al residues
            for residue in chain:
                res_name = residue.get_resname()
                if res_name not in [*three_to_one.keys()]:
                    print(res_name)
                    continue
                if is_aa(residue)!=True:
                    continue
                #Save residue
                model_seqs[chain.id]+=three_to_one[res_name]


    return model_seqs



def write_fasta(seq, seq_id, outname):
    """Write fasta
    """
    with open(outname, 'w') as file:
        file.write('>'+seq_id+'\n')
        file.write(seq)

##################MAIN#######################

#Parse args
args = parser.parse_args()
#Data
structure_file = args.structure_file[0]
outdir = args.outdir[0]
#Get seqs
model_seqs = read_pdb(structure_file)
pdbid = structure_file.split('/')[-1].split('.')[0]
#Write
for chain in model_seqs:
    write_fasta(model_seqs[chain], pdbid+'_'+chain, outdir+pdbid+'_'+chain+'.fasta')
