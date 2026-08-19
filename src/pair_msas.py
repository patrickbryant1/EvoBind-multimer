import argparse
import sys
import os
import numpy as np
import time
import pdb

parser = argparse.ArgumentParser(description = '''A simple script for matching sequences from two different MSAs based
                                                on the OX identifiers resulting from an hhblits search..''')

parser.add_argument('--a3m1', nargs=1, type= str, default=sys.stdin, help = 'Path to msa1 in a3m format.')
parser.add_argument('--a3m2', nargs=1, type= str, default=sys.stdin, help = 'Path to msa2 in a3m format.')
parser.add_argument('--max_gap_fraction', nargs=1, type=float, default=sys.stdin, help = 'The maximal gap fraction allowed in each sequence (default = 0.9).')
parser.add_argument('--outname', nargs=1, type= str, default=sys.stdin, help = 'Path to output filename')


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
                        try:
                            species.append(int(OX.split(' ')[0]))
                        except:
                            species.append(0)
                    else:
                        species.append(0)
                elif 'TaxID=' in line:
                    OX= line.split('TaxID=')[1]
                    if len(OX)>0:
                        try:
                            species.append(int(OX.split(' ')[0]))
                        except:
                            species.append(0)
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

def match_top_ox(ox1, ox2, msa1, msa2):
    '''Select the top ox match (first match) in each MSA and merge the sequences to a final MSA file
    in a3m format
    - number of possible combinations
    - median and std in hits per species
    '''
    #Don't remove the zeros (no OX), then the query sequences (first line)
    #will be removed
    matching_ox = np.intersect1d(ox1,ox2)

    ind1 = [] #Index to select from the individual MSAs
    ind2 = []
    ncombos = []
    #Go through all matching and select the first (top) hit
    for ox in matching_ox:
        ind1.append(min(np.argwhere(ox1==ox)[:,0]))
        ind2.append(min(np.argwhere(ox2==ox)[:,0]))

        ncombos.append(np.argwhere(ox1==ox).shape[0]*np.argwhere(ox2==ox).shape[0])

    #Select from MSAs and merge
    merged = np.concatenate((msa1[ind1], msa2[ind2]),axis=1)

    return merged, np.sum(ncombos), np.median(ncombos), np.std(ncombos)

def write_a3m(merged_msa, outfile):
    '''Write a3m MSA'''
    backmap = { 1:'A', 2:'C', 3:'D', 4:'E', 5:'F',6:'G' ,7:'H',
               8:'I', 9:'K', 10:'L', 11:'M', 12:'N', 13:'P',14:'Q',
               15:'R', 16:'S', 17:'T', 18:'V', 19:'W', 20:'Y', 21:'-'} #Here all unusual AAs and gaps are set to the same char (same in the GaussDCA script)

    with open(outfile,'w') as file:
        for i in range(len(merged_msa)):
            file.write('>'+str(i)+'\n')
            file.write(''.join([backmap[ch] for ch in merged_msa[i]])+'\n')

    return None
#################MAIN####################

#Parse args
args = parser.parse_args()
#Data
a3m1 = args.a3m1[0]
a3m2 = args.a3m2[0]
max_gap_fraction = args.max_gap_fraction[0]
outname = args.outname[0]

t1 = time.time()
#MSA1
msa1, ox1 = read_a3m(a3m1, max_gap_fraction)
#MSA2
msa2, ox2 = read_a3m(a3m2, max_gap_fraction)

#Get some statistics for msa1
nseqs_total1, l1 = msa1.shape
nunique_ox1 = np.unique(ox1).shape[0]-1 #Subtract 1 (no species)
#Get some statistics for msa2
nseqs_total2, l2 = msa2.shape
nunique_ox2 = np.unique(ox2).shape[0]-1 #Subtract 1 (no species)

#Match
merged_msa, ncombos_total, median_combos, std_combos = match_top_ox(ox1, ox2, msa1, msa2)

#Write the new a3m
write_a3m(merged_msa, outname)
#Print some statistics
t2 = time.time()
print('Saved paired MSA to', outname)
print('Matching and writing took', np.round(t2-t1,2), 'seconds.\nStatistics')
#print(id1+'\t'+id2)
print('Total seqs\t'+str(nseqs_total1)+'\t'+str(nseqs_total2))
print('Lengths\t'+str(l1)+'\t'+str(l2))
print('Unique organisms\t'+str(nunique_ox1)+'\t'+str(nunique_ox2))
print('-------------------------------------------------------------------------------')
print('Total combos\t'+str(ncombos_total))
print('Median combos\t'+str(median_combos))
print('Std combos\t'+str(std_combos))
print('Matched combos\t'+str(merged_msa.shape[0]))
print('Length of merged MSA\t'+str(merged_msa.shape[1]))
