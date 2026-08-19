import argparse
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import matplotlib.ticker as ticker
import glob
from collections import Counter
from scipy.stats import pearsonr, spearmanr
import scipy
import pdb

parser = argparse.ArgumentParser(description = '''Visualise anc calculate matrix.''')

parser.add_argument('--vhl_kras_metric_dir', nargs=1, type= str, default=sys.stdin, help = 'Path to dir with csv scores from design.')
parser.add_argument('--vhl_kras_experimental_results', nargs=1, type= str, default=sys.stdin, help = 'Path to csv with experimental results.')
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

def metric_vis(evobind_metrics, uid, outdir):
    """Plot metrics
    """

    def set_axis_style(ax, labels):
        ax.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
        ax.set_xlim(0.25, len(labels) + 0.75)

    print('Min length designed for '+uid+':', evobind_metrics.length.min())
    print('Max length designed for '+uid+':', evobind_metrics.length.max())

    #Length and loss
    top10 = []
    fig, ax = plt.subplots(figsize=(15/2.54, 9/2.54))
    for length in np.sort(evobind_metrics.length.unique()):
        sel = evobind_metrics[evobind_metrics.length==length]
        top10.append(sel.sort_values(by='loss').reset_index().loc[:int(len(sel)*0.1)-1].loss.values)

    plt.violinplot(top10,
                  showmeans=False,
                  showmedians=True)

    labels = np.sort(evobind_metrics.length.unique())
    ax.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
    tick_spacing = 1
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    plt.title('Loss distribution and length')
    plt.ylabel('Loss')
    plt.xlabel('Length')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outdir+uid+'_sequence_length_loss_distr.png', dpi=300, format='png')
    plt.close()


    #Optimisation curve per length
    fig, ax = plt.subplots(figsize=(15/2.54, 9/2.54))
    for length in np.sort(evobind_metrics.length.unique()):
        sel = evobind_metrics[evobind_metrics.length==length]
        min_sel = []
        for run in range(1,6):
            sel_run = sel[sel.run==run]
            min_sel.append(np.minimum.accumulate(sel_run.loss.values, axis=0))
        min_sel = np.array(min_sel)
        #Vis
        means = np.mean(min_sel,axis=0)
        stds = np.std(min_sel,axis=0)
        plt.plot(range(min_sel.shape[1]), means, label=str(length))
        #Fill std
        plt.fill_between(range(min_sel.shape[1]), means-stds, means+stds, alpha=0.15)

    plt.title('Optimisation curves')
    plt.yscale('log')
    plt.ylabel('Loss')
    plt.xlabel('Iteration')
    plt.legend()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outdir+uid+'_opt.png', dpi=300, format='png')
    plt.close()



    #Length and plDDT
    top10 = []
    fig, ax = plt.subplots(figsize=(15/2.54, 9/2.54))
    for length in np.sort(evobind_metrics.length.unique()):
        sel = evobind_metrics[evobind_metrics.length==length]
        top10.append(sel.sort_values(by='loss').reset_index().loc[:int(len(sel)*0.1)-1].plddt.values)
    plt.violinplot(top10,
                  showmeans=False,
                  showmedians=True)

    labels = np.sort(evobind_metrics.length.unique())
    ax.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
    tick_spacing = 1
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    plt.title('plDDT distribution and length')
    plt.ylabel('plDDT')
    plt.xlabel('Length')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outdir+uid+'sequence_length_plddt_distr.png', dpi=300, format='png')
    plt.close()

    #Top 10% overall coloured by plDDT
    fig, ax = plt.subplots(figsize=(12/2.54, 9/2.54))
    sel = evobind_metrics.sort_values(by='loss').reset_index().loc[:int(len(evobind_metrics)*0.1)-1]
    #Filter on plDDT>85
    sel = sel[sel.plddt>=85]
    print('Top 10% lowest loss with plDDT>=85:', len(sel))
    plt.scatter(sel.plddt, sel.loss, c=sel.length, s=1, marker='+', alpha=0.5)
    plt.colorbar(label='Length')
    plt.title('Top 10% coloured by plDDT')
    plt.ylabel('Loss')
    plt.xlabel('plDDT')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outdir+uid+'length_vs_loss_plddt_c.png', dpi=300, format='png')
    plt.close()


def experimental_vis(vhl_kras_experimental_results, outdir):
    """Visualise the experimental results
    """


    pos_ctrl_conc = [0, 0.00000953, 0.0000381, 0.00015258, 0.000610, 0.00244, 0.009765625, 0.0390625, 0.15625, 0.625, 2.5, 10]
    pep_conc = [0, 9.766, 19.531, 39.0625, 78.125, 156.25, 312.5, 625, 1250, 2500, 5000, 10000]


    #Subtract background
    for replicate in ['1', '2']:
        for col in ['LC-2_', 'peptide_743_', 'peptide_587_', 'peptide_587_']:
            vhl_kras_experimental_results[col+replicate] = vhl_kras_experimental_results[col+replicate]-vhl_kras_experimental_results['background_'+replicate]

    #Vis mean and std
    pretty_names =  {'peptide_743_':'RPGDPVCSWW'}
    IC50 = {'peptide_743_':'2082 μM'}
    colors = {'peptide_743_':'magenta'}
    concentrations = {'peptide_743_':pep_conc}
    fig, ax = plt.subplots(figsize=(9/2.54, 9/2.54))
    for col in pretty_names:
        cat_vals = np.concatenate([np.array([vhl_kras_experimental_results[col+'1']]), np.array([vhl_kras_experimental_results[col+'2']])],axis=0)
        means = np.mean(cat_vals,axis=0)
        stds = np.std(cat_vals,axis=0)
        plt.plot(concentrations[col], means, label='IC50='+IC50[col], color=colors[col])
        plt.fill_between(concentrations[col], means-stds, means+stds, alpha=0.25, color=colors[col])

    plt.title('RPGDPVCSWW')
    plt.ylabel('RLU')
    plt.xlabel('μM')
    plt.xscale('log')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir+'vhl_kras_deg_RPGDPVCSWW.png', dpi=300, format='png')
    plt.close()

    pretty_names =  {'peptide_587_':'APGDPVCSWW'}
    IC50 = {'peptide_587_':'1568 μM'}
    colors = {'peptide_587_':'magenta'}
    concentrations = {'peptide_587_':pep_conc}
    fig, ax = plt.subplots(figsize=(9/2.54, 9/2.54))
    for col in pretty_names:
        cat_vals = np.concatenate([np.array([vhl_kras_experimental_results[col+'1']]), np.array([vhl_kras_experimental_results[col+'2']])],axis=0)
        means = np.mean(cat_vals,axis=0)
        stds = np.std(cat_vals,axis=0)
        plt.plot(concentrations[col], means, label='IC50='+IC50[col], color=colors[col])
        plt.fill_between(concentrations[col], means-stds, means+stds, alpha=0.25, color=colors[col])

    plt.title('APGDPVCSWW')
    plt.ylabel('RLU')
    plt.xlabel('μM')
    plt.xscale('log')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir+'vhl_kras_deg_APGDPVCSWW.png', dpi=300, format='png')
    plt.close()

    pretty_names =  {'LC-2_':'+ctrl', 'peptide_587_':'APGDPVCSWW', 'peptide_743_':'RPGDPVCSWW'}
    colors = {'LC-2_':'grey', 'peptide_587_':'tab:green', 'peptide_743_':'tab:blue'}
    concentrations = {'LC-2_':pos_ctrl_conc, 'peptide_587_':pep_conc, 'peptide_743_':pep_conc}
    fig, ax = plt.subplots(figsize=(9/2.54, 9/2.54))
    for col in pretty_names:
        cat_vals = np.concatenate([np.array([vhl_kras_experimental_results[col+'1']]), np.array([vhl_kras_experimental_results[col+'2']])],axis=0)
        means = np.mean(cat_vals,axis=0)
        stds = np.std(cat_vals,axis=0)
        plt.plot(concentrations[col], means, label=pretty_names[col], color=colors[col])
        plt.fill_between(concentrations[col], means-stds, means+stds, alpha=0.25, color=colors[col])

    plt.title('Degradation curves')
    plt.ylabel('RLU')
    plt.xlabel('μM')
    plt.xscale('log')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir+'vhl_kras_deg_all.png', dpi=300, format='png')
    plt.close()

    pdb.set_trace()


##################MAIN#######################

#Parse args
args = parser.parse_args()
#Data
vhl_kras_metric_dir = args.vhl_kras_metric_dir[0]
vhl_brd4_metric_dir = args.vhl_brd4_metric_dir[0]
vhl_kras_experimental_results = pd.read_csv(args.vhl_kras_experimental_results[0])
outdir = args.outdir[0]

# vhl_kras_metrics = merge_metrics(vhl_kras_metric_dir)
# vhl_brd4_metrics = merge_metrics(vhl_brd4_metric_dir)
#
# metric_vis(vhl_kras_metrics, 'vhl_kras', outdir)
# metric_vis(vhl_brd4_metrics, 'vhl_brd4', outdir)
experimental_vis(vhl_kras_experimental_results, outdir)

pdb.set_trace()
