#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.cbook import boxplot_stats
from matplotlib.backends.backend_pdf import PdfPages


def get_args():
    parser = argparse.ArgumentParser(description='Evaluate differntial binding')
    parser.add_argument('file1', metavar='CONDITION1.tsv', type=str,
                    help='TSV file from pipeline run with first condition')
    parser.add_argument('file2', metavar='CONDITION2.tsv', type=str,
                    help='TSV file from pipeline run with seccond condition')
    parser.add_argument('--BD-comp',default=None, metavar="BINDetect_results.txt", type=str,help="get a compairison Plot to BINDetect results" )
    parser.add_argument('--out', metavar="FILENAME_prefix",help="Prefix of how the output files should be named.", required=True)
    parser.add_argument('--get-outliers', action='store_true', help="Set this flag for an additional file containing the names of the motifs in the lowest and highest quartile.")
    parser.add_argument('--cond_names', nargs=2, required=True, help="Names of the two conditions that are compared. If --BD-comp flag ist set, pleas make sure that the same condition names were used in BINDetect")
    args = parser.parse_args()
    return args


def plots_and_stats(df_comp,df_all,cond_names,args):
    # create plots:
    # Plotting
    sns.set_theme(style="whitegrid")
    pdf = PdfPages(f"{args.out}.pdf")
    if args.BD_comp:
        scat = sns.regplot(data=df_all, x="DuxPos_DuxNeg_change",y="ratio_change")
        scat.set_xlabel('BINDetect Change')
        scat.set_ylabel('Pipeline ratio change')
        props = dict(boxstyle='square', facecolor='lightsteelblue', alpha=0.5)
        scat.text(0,1.7,f'r = {round(np.corrcoef(df_all[f"{cond_names[0]}_{cond_names[1]}_change"].to_numpy(),df_all["ratio_change"].to_numpy())[0,1],2)}',
                bbox=props, size=20)
        pdf.savefig(scat.get_figure())
        plt.clf()

    violin = sns.violinplot( x="ratio_change", data=df_comp)
    violin.set_xlabel('Change in ratio')
    pdf.savefig(violin.get_figure())
    plt.clf()

    bp = sns.boxplot(x="ratio_change", data=df_comp)
    bp.set_xlabel('Change in ratio')
    pdf.savefig(bp.get_figure())
    plt.clf()
    pdf.close()
    
    stats =  boxplot_stats(df_comp['ratio_change'])[0] if args.get_outliers else None

    return stats



def main():
    args = get_args()
    # read in files
    cond_names = args.cond_names
    file_name1, file_name2 = args.file1, args.file2 
    bindetect_file = args.BD_comp
    
    print('reading in and filtering files...')
    df1 = pd.read_csv(file_name1,sep='\t',usecols=['id','name','nsites_whole_genome','nsites_open_chromatin'])
    df2 = pd.read_csv(file_name2,sep='\t', usecols=['id','nsites_whole_genome','nsites_open_chromatin']) 

    # filter motifs of seccond run:
    df1 = df1[~df1.id.str.contains('^motif_\d+$')]
    df2 = df2[~df2.id.str.contains('^motif_\d+$')]

    print('calculating percentage of open chromatin...')
    # calculate percentage open
    df1["p_open_c"] = df1["nsites_open_chromatin"]*100/df1["nsites_whole_genome"]
    df2["p_open_c"] = df2["nsites_open_chromatin"]*100/df2["nsites_whole_genome"]

    df_merged = df1.join(df2.set_index('id'),how='inner',lsuffix=f"_{cond_names[0]}", rsuffix=f"_{cond_names[1]}", on='id')
    df_merged['ratio_change'] = df_merged[f'p_open_c_{cond_names[0]}'] - df_merged[f'p_open_c_{cond_names[1]}']

    # add Bind detect if applicable
    df_all = None
    if bindetect_file:
        print('loading BINDetect file...')
        df_BD = pd.read_csv(bindetect_file, sep='\t', usecols=["motif_id",f"{cond_names[0]}_{cond_names[1]}_change"])
        df_all = df_merged.join(df_BD.set_index('motif_id'), how='inner',lsuffix='_in',rsuffix='_BD', on='id')

    # Plotting and extracting stats
    print('plotting...')
    stats = plots_and_stats(df_merged, df_all, cond_names, args)

    # interesting motifs
    if stats:
        print('extracting interesting motifs...')
        with open(f"{args.out}_outlier_motifs.txt",'w') as out_file:
            out_file.write('> Lower Quartile:\n')
            df_motf = df_all[df_merged['ratio_change'] < stats["q1"]]
            for name in df_motf["name"]:
                out_file.write(f'{name}\n')

            out_file.write('> Upper Quartile:\n')
            df_motf = df_all[df_merged['ratio_change'] > stats["q3"]]
            for name in df_motf["name"]:
                out_file.write(f'{name}\n')
    print('Done!')

if __name__ == "__main__":
    main()
