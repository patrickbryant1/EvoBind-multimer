import argparse
import sys
import os
import numpy as np
import pandas as pd
import glob
import pdb

parser = argparse.ArgumentParser(description = '''Fuse the MSAs for interacting chains by writing gaps where the other chain should be.''')

parser.add_argument('--a3m1', nargs=1, type= str, default=sys.stdin, help = 'Path to a3m file 1.')
parser.add_argument('--a3m2', nargs=1, type= str, default=sys.stdin, help = 'Path to a3m file 2.')
parser.add_argument('--max_gap_fraction', nargs=1, type=float, default=sys.stdin, help = 'The maximal gap fraction allowed in each sequence (default = 0.9).')
parser.add_argument('--outname', nargs=1, type= str, default=sys.stdin, help = 'Path to file to write to.')

def read_a3m(infile,max_gap_fraction=0.9):
    '''Read a3m MSA'''
    mapping = {'-': 21, 'A': 1, 'B': 21, 'C': 2, 'D': 3, 'E': 4, 'F': 5,
             'G': 6,'H': 7, 'I': 8, 'K': 9, 'L': 10, 'M': 11,'N': 12,
             'O': 21, 'P': 13,'Q': 14, 'R': 15, 'S': 16, 'T': 17,
             'V': 18, 'W': 19, 'Y': 20,'U': 21, 'Z': 21, 'X': 21, 'J': 21}

    parsed = []#Save extracted msa
    species = []
    seqlen = 0
    lc = 0
    with open(infile, 'r') as file:
        for line in file:
            line = line.rstrip()

            if line.startswith('>'): #OX=OrganismIdentifier
                if 'OX=' in line:
                    OX= line.split('OX=')[1]
                    if len(OX)>0:
                        species.append(int(OX.split(' ')[0]))
                    else:
                        species.append(0)
                else:
                    species.append(0)
                continue
            line = line.rstrip()
            gap_fraction = line.count('-') / float(len(line))
            if gap_fraction <= max_gap_fraction:#Only use the lines with less than 90 % gaps
                parsed.append([mapping.get(ch, 22) for ch in line if not ch.islower()])
            else:
                if len(species)>1:
                    species = species[:-1] #Remove the previously stored species
                    continue
            #Check that the lengths match
            if len(parsed[-1])!=seqlen and lc>=1:
                parsed = parsed[:-1]
                species = species[:-1]
                continue
            seqlen = len(parsed[-1])
            lc+=1


    return np.array(parsed, dtype=np.int8, order='F'), np.array(species)

def write_a3m(fused, outfile):
    '''Write a3m MSA'''
    backmap = { 1:'A', 2:'C', 3:'D', 4:'E', 5:'F',6:'G' ,7:'H',
               8:'I', 9:'K', 10:'L', 11:'M', 12:'N', 13:'P',14:'Q',
               15:'R', 16:'S', 17:'T', 18:'V', 19:'W', 20:'Y', 21:'-'} #Here all unusual AAs and gaps are set to the same char (same in the GaussDCA script)

    with open(outfile,'w') as file:
        for i in range(len(fused)):
            file.write('>'+str(i)+'\n')
            file.write(''.join([backmap[ch] for ch in fused[i]])+'\n')

    return None


#################MAIN####################

#Parse args
args = parser.parse_args()
max_gap_fraction = args.max_gap_fraction[0]
#Data
a3m1, species1 = read_a3m(args.a3m1[0], max_gap_fraction)
a3m2, species2 = read_a3m(args.a3m2[0], max_gap_fraction)
outname = args.outname[0]

#Construct entire a3m matrix
fused = np.zeros((a3m1.shape[0]+a3m2.shape[0],a3m1.shape[1]+a3m2.shape[1]))
fused[:]=21 #Assign gaps
#Assign a3m1
fused[:a3m1.shape[0],:a3m1.shape[1]]=a3m1
#Assign a3m2
fused[a3m1.shape[0]:,a3m1.shape[1]:]=a3m2
#Write the fused MSA
write_a3m(fused, outname)
print('Saved paired MSA to', outname)
