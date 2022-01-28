#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.cbook import boxplot_stats
from matplotlib.backends.backend_pdf import PdfPages
# read filenames
parser = argparse.ArgumentParser(description='Evaluate differntial binding')
parser.add_argument('file1')
parser.add_argument('file2')
parser.add_argument('--BD-comp')
parser.add_argument('--out')
parser.add_argument('--get-outliers')
parser.add_argument('--cond_names')
args = parser.parse_args()

# reas in files
cond_names = ["cond1","cond2"]
file_name1, file_name2 = "motif_ranks_Pos_sub.tsv","motif_ranks_Neg_sub.tsv" #sys.argv[1:3]
bindetect_file = "bindetect_results_TOBIAS.txt"

df1 = pd.read_csv(file_name1,sep='\t',usecols=['id','name','nsites_whole_genome','nsites_open_chromatin']) #enrichment_score
df2 = pd.read_csv(file_name2,sep='\t', usecols=['id','nsites_whole_genome','nsites_open_chromatin']) # enrichment_score
df_BD = pd.read_csv(bindetect_file, sep='\t', usecols=["motif_id",f"{cond_names[0]}_{cond_names[1]}_change"])

# filter motifs of seccond run:
df1 = df1[~df1.id.str.contains('^motif_\d+$')]
df2 = df2[~df2.id.str.contains('^motif_\d+$')]


# calculate percentage open
df1["p_open_c"] = df1["nsites_open_chromatin"]*100/df1["nsites_whole_genome"]
df2["p_open_c"] = df2["nsites_open_chromatin"]*100/df2["nsites_whole_genome"]

df_merged = df1.join(df2.set_index('id'),how='inner',lsuffix=cond_names[0], rsuffix=cond_names[1], on='id')
df_merged['ratio_change'] = df_merged[f'p_open_c_{cond_names[0]}'] - df_merged[f'p_open_c_{cond_names[1]}']

df_all = df_merged.join(df_BD.set_index('motif_id'), how='inner',lsuffix='_in',rsuffix='_BD', on='id')

# Plotting
sns.set_theme(style="whitegrid")
scat = sns.regplot(data=df_all, x="DuxPos_DuxNeg_change",y="ratio_change")
scat.set_xlabel('BINDetect Change')
scat.set_ylabel('Pipeline ratio change')
props = dict(boxstyle='square', facecolor='lightsteelblue', alpha=0.5)
scat.text(0,1.7,f'r = {round(np.corrcoef(df_all[f"{cond_names[0]}_{cond_names[1]}_change"].to_numpy(),df_all["ratio_change"].to_numpy())[0,1],2)}',
          bbox=props, size=20)
plt.show()

violin = sns.violinplot( x="ratio_change", data=df_merged)
violin.set_xlabel('Change in ratio')
plt.show()

bp = sns.boxplot(x="ratio_change", data=df_merged)
bp.set_xlabel('Change in ratio')
plt.show()



pdf = PdfPages("test_figs.pdf")
pdf.savefig(bp.get_figure())
pdf.savefig(violin.get_figure())
pdf.savefig(scat.get_figure())
pdf.close()

stats = boxplot_stats(df_merged['ratio_change'])[0]

# interesting motifs
with open('motifs_outer_quartiles.txt','w') as out_file:
    out_file.write('> Lower Quartile:\n')
    df_motf = df_all[df_merged['ratio_change'] < stats["q1"]]
    for name in df_motf["name"]:
        out_file.write(f'{name}\n')

    out_file.write('> Upper Quartile:\n')
    df_motf = df_all[df_merged['ratio_change'] > stats["q3"]]
    for name in df_motf["name"]:
        out_file.write(f'{name}\n')
