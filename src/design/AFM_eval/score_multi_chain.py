import os
import pickle
import sys
import time
import argparse
import pandas as pd
import numpy as np
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SVDSuperimposer import SVDSuperimposer

parser = argparse.ArgumentParser(description="Score AFM vs EvoBind.")

parser.add_argument('--afm_dir', type=str, required=True, help='Path to location of AFM preds.')
parser.add_argument('--evo_dir', type=str, required=True, help='Path to location of EvoBind preds.')
parser.add_argument('--score_csv', type=str, required=True, help='Path to csv with sequences and ids to be scores.')
parser.add_argument('--outdir', type=str, required=True, help='Path to output directory. Include / in end')

def read_pdb(pdbname):
    '''Read PDB
    '''

    f=open(pdbname,'rt')

    if '.pdb' in pdbname:
        parser = PDBParser()
        struc = parser.get_structure('', f)
    else:
        parser = MMCIFParser()
        struc = parser.get_structure('',f)

    #Save
    model_coords = {}
    model_3seq = {}
    model_resnos = {}
    model_atoms = {}
    model_bfactors = {}


    for model in struc:
        for chain in model:
            #Save
            model_coords[chain.id]=[]
            model_3seq[chain.id]=[]
            model_resnos[chain.id]=[]
            model_atoms[chain.id]=[]
            model_bfactors[chain.id]=[]

            #Go through al residues
            for residue in chain:
                res_name = residue.get_resname()
                if is_aa(residue)!=True:
                    continue
                for atom in residue:
                    atom_id = atom.get_id()
                    atm_name = atom.get_name()
                    #Save
                    model_coords[chain.id].append(atom.get_coord())
                    model_3seq[chain.id].append(res_name)
                    model_resnos[chain.id].append(residue.get_id()[1])
                    model_atoms[chain.id].append(atom_id)
                    model_bfactors[chain.id].append(atom.bfactor)



    return model_coords, model_3seq, model_resnos, model_atoms, model_bfactors

def score_by_target_chain(target_ci, peptide_ci, evo_chains, afm_chains,
        evo_coords, evo_3seq, evo_resnos, evo_atoms, evo_bfactors,
        afm_coords, afm_3seq, afm_resnos, afm_atoms, afm_bfactors):
    
    # Define coords
    evo_rec_coords, evo_pep_coords = np.array(evo_coords[evo_chains[target_ci]]), np.array(evo_coords[evo_chains[peptide_ci]])
    afm_rec_coords, afm_pep_coords = np.array(afm_coords[afm_chains[target_ci]]), np.array(afm_coords[afm_chains[peptide_ci]])

    # Align CAs
    sup = SVDSuperimposer()
    evo_rec_CAs = evo_rec_coords[np.argwhere(np.array(evo_atoms[evo_chains[target_ci]])=='CA')[:,0]]
    afm_rec_CAs = afm_rec_coords[np.argwhere(np.array(afm_atoms[afm_chains[target_ci]])=='CA')[:,0]]
    
    # Ensure the number of CAs match
    min_length = min(len(evo_rec_CAs), len(afm_rec_CAs))
    evo_rec_CAs = evo_rec_CAs[:min_length]
    afm_rec_CAs = afm_rec_CAs[:min_length]
    
    # Check if we have enough points to perform the superimposition
    if min_length < 3:
        print(f"Warning: Not enough CA atoms to perform superimposition. EvoBind: {len(evo_rec_CAs)}, AFM: {len(afm_rec_CAs)}")
        return None, None, None

    sup.set(evo_rec_CAs, afm_rec_CAs)
    sup.run()
    rot, tran = sup.get_rotran()
    
    # Rotate the peptide coords to match the centre of mass for the native comparison
    evo_pep_CAs = evo_pep_coords[np.argwhere(np.array(evo_atoms[evo_chains[peptide_ci]])=='CA')[:,0]]
    afm_pep_CAs = afm_pep_coords[np.argwhere(np.array(afm_atoms[afm_chains[peptide_ci]])=='CA')[:,0]]
    rotated_afm_pep_CAs = np.dot(afm_pep_CAs, rot) + tran
    
    # Peptide COM
    evo_pep_CM = np.mean(evo_pep_CAs, axis=0)
    afm_pep_CM = np.mean(rotated_afm_pep_CAs, axis=0)
    delta_CM = np.linalg.norm(evo_pep_CM - afm_pep_CM)

    # Get EvoBind if pos - CBs<8Å
    evo_rec_CBs = evo_rec_coords[np.argwhere(np.array(evo_atoms[evo_chains[target_ci]])=='CB')[:,0]]
    evo_pep_CBs = evo_pep_coords[np.argwhere(np.array(evo_atoms[evo_chains[peptide_ci]])=='CB')[:,0]]
    mat = np.vstack((evo_rec_CBs, evo_pep_CBs))
    evo_dmat = np.linalg.norm(mat[:, np.newaxis, :] - mat[np.newaxis, :, :], axis=2)
    
    # Get interface
    l1 = len(evo_rec_CBs)
    contact_dists = evo_dmat[:l1, l1:]
    rec_if_pos = np.unique(np.argwhere(contact_dists < 8)[:, 0])
    if_resno_pos = []
    for ri in rec_if_pos:
        if_resno_pos.extend(np.argwhere(np.array(evo_resnos[evo_chains[target_ci]]) == ri+1)[:, 0])

    # Get AFM interface dists
    if len(if_resno_pos) > 0:
        afm_if_coords = np.vstack((afm_rec_coords[if_resno_pos], afm_pep_coords))
        afm_dmat = np.linalg.norm(afm_if_coords[:, np.newaxis, :] - afm_if_coords[np.newaxis, :, :], axis=2)
        afm_if = afm_dmat[:len(if_resno_pos), len(if_resno_pos):]
        
        afm_rec_avg_if_dist = np.mean(np.min(afm_if, axis=1))
        afm_pep_avg_if_dist = np.mean(np.min(afm_if, axis=0))
    else:
        afm_rec_avg_if_dist = 20
        afm_pep_avg_if_dist = 20
    
    # Get the peptide plDDT
    afm_pep_plddt = np.mean(afm_bfactors[afm_chains[peptide_ci]])
    
    # Calc AFM-EvoBind loss
    loss = (afm_rec_avg_if_dist + afm_pep_avg_if_dist) / 2 * 1 / afm_pep_plddt * delta_CM
    
    return loss, afm_pep_plddt, delta_CM

def score_structures(afm_dir, evo_dir, score_df, outdir):
    loss_df = {'iter': [], 'sequence': [], 'loss': [], 'plddt': [], 'COM': []}

    for _, row in score_df.iterrows():
        if row.iteration == 'init':
            continue

        try:
            evo_coords, evo_3seq, evo_resnos, evo_atoms, evo_bfactors = read_pdb(os.path.join(evo_dir, f'unrelaxed_{row.iteration}.pdb'))
        except Exception as e:
            print(f'Could not read file: {os.path.join(evo_dir, f"unrelaxed_{row.iteration}.pdb")}')
            print(f'Error: {str(e)}')
            continue

        try:
            afm_coords, afm_3seq, afm_resnos, afm_atoms, afm_bfactors = read_pdb(os.path.join(afm_dir, f'{row.sequence}.pdb'))
        except Exception as e:
            print(f'Could not read file: {os.path.join(afm_dir, f"{row.sequence}.pdb")}')
            print(f'Error: {str(e)}')
            continue

        evo_chains = list(evo_coords.keys())
        afm_chains = list(afm_coords.keys())

        result = score_by_target_chain(
            0, 2, evo_chains, afm_chains,
            evo_coords, evo_3seq, evo_resnos, evo_atoms, evo_bfactors,
            afm_coords, afm_3seq, afm_resnos, afm_atoms, afm_bfactors
        )
        
        if result is None:
            print(f"Skipping iteration {row.iteration} due to insufficient data for superimposition.")
            continue
        
        loss, afm_pep_plddt, delta_CM = result

        loss_df['iter'].append(str(row.iteration))
        loss_df['sequence'].append(row.sequence)
        loss_df['loss'].append(loss)
        loss_df['plddt'].append(afm_pep_plddt)
        loss_df['COM'].append(delta_CM)

    loss_df = pd.DataFrame(loss_df)
    output_file = os.path.join(outdir, 'afm_evo_metrics.csv')
    loss_df.to_csv(output_file, index=None)
    print(f'Saved losses to {output_file}')

if __name__ == "__main__":
    args = parser.parse_args()
    score_df = pd.read_csv(args.score_csv)
    score_structures(args.afm_dir, args.evo_dir, score_df, args.outdir)