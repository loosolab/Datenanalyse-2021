#!/usr/bin/env python3

import pandas as pd
import plotly.express as px
import argparse
import re

parser = argparse.ArgumentParser(description='Evaluate differntial binding')
parser.add_argument('--motifs', metavar='MOTIF_CLUSER.yml', help='Output of the motif clustering with TOBIAS', required=True)
parser.add_argument('--out', metavar="FILENAME_prefix",help="Prefix of how the output files should be named.", required=True)
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
        tissue, cell_type = cur_motif_name.split('_')
        tissues.append(tissue)
        cell_types.append(re.sub('\d+$', '', cell_type))

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

### Plotting###
## create Dotplot
fig = px.scatter(df_motifs, x="Cell_type", y="Tissue", color="Cluster", hover_data=['motif_name'], title="Motif Similarity relations./run")
fig.update_traces(marker=dict(size=15),
                  selector=dict(mode='markers'))
fig.show()
# TODO filepath
fig.write_html(f"/mnt/workspace_stud/allstud/wp5/runs/{args.out}_Dotplot.html")

## create count bar plot
# add count to Data Frame
df_motifs["count"] = pd.Series([1]*df_motifs.shape[0])
# plot with plotly
fig = px.bar(df_motifs, x="Cluster", y="count", color="Tissue", title="Motif Counts for all new Motifs")
fig.update_layout(barmode='stack', xaxis={'categoryorder':'total descending'})
fig.update_traces(marker_line_width=0)
fig.show()
# TODO filepath
fig.write_html(f"/mnt/workspace_stud/allstud/wp5/runs/{args.out}_bar.html")

