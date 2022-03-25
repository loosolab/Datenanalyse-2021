#!/usr/bin/env python3

import pandas as pd
import csv
import os

## prepare Dataframe 
# initialize helping lists
cur_cluster = ''
motif_names = []
clusters = []
tissues = []
cell_types = []

# iterate over files
files = os.listdir("Desktop/myFolder")
myfile = 'filename.txt'

for filename in files:
    if filename == myfile:
        continue


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
