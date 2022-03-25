#!/usr/bin/env python3

import pandas as pd
import plotly.express as px
import argparse
import csv

parser = argparse.ArgumentParser(description='Evaluate differntial binding')
parser.add_argument('--runs-dir', metavar="PATH/TO/RUNS", help='Path to directory where the motif discovery runs are stored', required=True )
parser.add_argument('--motifs', metavar='MOTIF_CLUSTER.yml', help='Output of the motif clustering with TOBIAS', required=True)
parser.add_argument('--out', metavar="FILENAME_prefix",help="Prefix of how the output files should be named.", required=True)
parser.add_argument('--annotation_dir', metavar="PATH/TO/WP2", help='Path to where the preproceesed original data lies. Should contain an annotation.txt for all tissues with the assignment of cluster to tissue.')
args = parser.parse_args()

## prepare Dataframe 
# initialize helping lists
cur_cluster = ''
motif_names = []
clusters = []
tissues = []
cell_types = []

# read file
with open(args.motifs) as file:
    for line in file:
        if "Cluster" in line:
            cur_cluster = line.strip()[:-1]
            continue
        cur_motif_name = line.strip().split(' ')[1]
        # drop known motifs
        if '.' in cur_motif_name:
            continue
        clusters.append(cur_cluster)
        motif_names.append(cur_motif_name)
        tissue, cell_type, _ = cur_motif_name.split('-')
        tissues.append(tissue)
        cell_types.append(cell_type)


# parse to dataframe
df_motifs = pd.DataFrame({
    'motif_name': motif_names,
    'Cluster': clusters,
    'Tissue': tissues,
    'Cell_type': cell_types,
})

# remove unneeded helping lists
del clusters
del tissues
del cell_types
del motif_names

# if available transalte cluster numbers to names
if args.annotation_dir :
    anno_dict = dict()
    for tissue in pd.unique(df_motifs['Tissue']):
        FILE=f"{args.annotation_dir}/{tissue}/annotation/annotation.txt"
        t_dict = dict()
        with open(FILE) as file:
            for line in csv.reader(file,delimiter='\t'):
                t_dict[f'cluster{line[0]}'] = line[1]
        anno_dict[tissue] = t_dict
    
    df_motifs['Cell_type'] = df_motifs.apply(lambda x: anno_dict[x['Tissue'][x['Cell_type']]]  ,axis=1) # TODO test
    
# write to csv
df_motifs.to_csv(f"{args.runs_dir}/{args.out}_clustering.csv")



### Plotting###
## create Dotplot
fig = px.scatter(df_motifs, x="Cell_type", y="Tissue", color="Cluster", hover_data=['motif_name'], title="Motif Similarity relations")
fig.update_traces(marker=dict(size=15),
                  selector=dict(mode='markers'))

# save
fig.write_html(f"{args.runs_dir}/{args.out}_Dotplot.html")

## create count bar plot
# add count to Data Frame
df_motifs["count"] = pd.Series([1]*df_motifs.shape[0])
# plot with plotly
fig = px.bar(df_motifs, x="Cluster", y="count", color="Tissue", title="Motif Counts for all new Motifs")
fig.update_layout(barmode='stack', xaxis={'categoryorder':'total descending'})
fig.update_traces(marker_line_width=0)

# save
fig.write_html(f"{args.runs_dir}/{args.out}_bar.html")

