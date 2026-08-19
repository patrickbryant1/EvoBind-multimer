import argparse
import sys
import os
import numpy as np
import pandas as pd
import glob
from collections import Counter
import pdb

parser = argparse.ArgumentParser(description = '''Select from the EvoBind optimisation.''')

parser.add_argument('--vhl_kras_metric_dir', nargs=1, type= str, default=sys.stdin, help = 'Path to dir with csv scores from design.')
parser.add_argument('--vhl_brd4_metric_dir', nargs=1, type= str, default=sys.stdin, help = 'Path to dir with csv scores from design.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')



##############FUNCTIONS##############
def merge_metrics(metric_dir):
    """Merge metrics from the opt
    """

    cat_df = []
    for name in glob.glob(metric_dir+'*.csv'):
        df = pd.read_csv(name)
        run = name.split('/')[-1].split('_')[1]
        length = name.split('/')[-1].split('_')[2].split('.')[0]
        df['run']=int(run)
        df['length']=int(length)
        cat_df.append(df)

    cat_df = pd.concat(cat_df)
    return cat_df

def select(evobind_metrics, uid):
    """
    Select
    1. Top 10% from loss
    2. plDDT >=85
    3. One per length
    """
    #Solubility
    #>25% charged residues (D, K, R, H and E) and <25% hydrophobic residues (W, L, I, F, M, V, Y).
    resi_counts = [Counter(x) for x in evobind_metrics.sequence]
    num_charged = [x['D']+x['K']+x['R']+x['H']+x['E'] for x in resi_counts]
    num_hydro = [x['W']+x['L']+x['I']+x['F']+x['M']+x['V']+x['Y'] for x in resi_counts]
    evobind_metrics['num_charged_resis']=num_charged
    evobind_metrics['num_hydro_resis']=num_hydro
    evobind_metrics['frac_charged']=np.array(num_charged)/evobind_metrics.length.values
    evobind_metrics['frac_hydro']=np.array(num_hydro)/evobind_metrics.length.values
    evobind_metrics['solubility_pass'] = 0
    evobind_metrics = evobind_metrics.reset_index()
    evobind_metrics.loc[evobind_metrics[(evobind_metrics.frac_charged>0.25)&(evobind_metrics.frac_hydro<0.25)].index,'solubility_pass']=1

    #Get top 10%
    sel = evobind_metrics.sort_values(by='loss').reset_index().loc[:int(len(evobind_metrics)*0.1)-1]
    #Filter on plDDT>85
    sel = sel[sel.plddt>=85]
    sel = sel[sel.columns[2:]]
    print('Selected top10% with plDDT>=85', len(sel))
    #Get top 100
    sel = sel.reset_index().loc[:99]
    print(Counter(sel.length))
    #Select unique lengths
    #sel = sel.loc[sel.length.drop_duplicates().index]
    print('Selected',len(sel), 'sequences for', uid)
    #Save
    sel[['run', 'iteration','sequence', 'length','if_dist_peptide_t1','if_dist_peptide_t2','plddt','loss','frac_charged','frac_hydro','solubility_pass']].to_csv(outdir+uid+'_top10_perc_plddt_over_85_u_length.csv', index=None)

##################MAIN#######################

#Parse args
args = parser.parse_args()
#Data
vhl_kras_metric_dir = args.vhl_kras_metric_dir[0]
vhl_brd4_metric_dir = args.vhl_brd4_metric_dir[0]
outdir = args.outdir[0]
#Merge
vhl_kras_metrics = merge_metrics(vhl_kras_metric_dir)
vhl_brd4_metrics = merge_metrics(vhl_brd4_metric_dir)
#Select
select(vhl_kras_metrics, 'vhl_kras')
select(vhl_brd4_metrics, 'vhl_brd4')
