#!/usr/bin/env python3
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter

# Initialize data lists for tissue, cell_type und motifs
data = pd.read_csv('/mnt/workspace_stud/allstud/wp5/runs/check_new_motifs_all.txt', sep='\t')

df = pd.DataFrame(data)

tissues = list(df.iloc[:, 0])
# cell_types = list(df.iloc[:, 1])  # reactivate to generate stats for known cell types
motifs = list(df.iloc[:, 2])

# check how often a motiv occurs
stats_all_motifs = Counter(motifs)
motif_counts = dict(stats_all_motifs)

# create plot for all tissues
# change the figure size
fig = plt.figure(figsize=(18, 14), dpi=300)

# plot the data as bar plot
plt.bar(motif_counts.keys(), motif_counts.values(), color='r')
plt.xticks(rotation=90)
plt.title("Motif Counts for all new Motifs (all tissues)")
plt.xlabel("New Motifs", labelpad=15)
plt.ylabel("Counts")

# show bar plot
plt.show()
    
# save plot
figure_name = "new_motif_counts_for_all_tissues_barplot.png"
path = "/mnt/workspace_stud/allstud/wp5/runs/"  # path to save the figure
fig.savefig(path + figure_name, dpi=fig.dpi)

# create plot per tissue
while len(tissues) != 0:

    # get tissue name and initialize count to get all cell types per tissue
    new_tissue = tissues[0]
    tissue_count = 1

    # check the quantity of different cell_types in a tissue
    for i in range(len(tissues)-1):
        if tissues[i] == tissues[i+1]:
            tissue_count += 1

    # initialize lists per tissue
    tissue_x = tissues[0:tissue_count]
    # cell_type_x = cell_types[0:tissue_count]  # reactivate to generate stats for known cell types
    motifs_x = motifs[0:tissue_count]

    # remove used values from data lists
    tissues = tissues[tissue_count:]
    #cell_types = cell_types[tissue_count:]  # reactivate to generate stats for known cell types
    motifs = motifs[tissue_count:]

    # check how often a motiv occurs
    stats_tissue_motifs = Counter(motifs_x)
    motif_tissue_counts = dict(stats_tissue_motifs)
    
    # change the figure size
    fig = plt.figure(figsize=(18, 14), dpi=300)

    # plot the data as bar plot
    plt.bar(motif_tissue_counts.keys(), motif_tissue_counts.values(), color='b')
    plt.xticks(rotation=90)
    plt.title("Motif Counts for all new Motifs (" + new_tissue + ")")
    plt.xlabel("New Motifs", labelpad=15)
    plt.ylabel("Counts")

    # show bar plot
    plt.show()
    
    # save plot
    figure_name = "new_motif_counts_for_" + new_tissue + "_barplot.png"
    path = "/mnt/workspace_stud/allstud/wp5/runs/" + new_tissue + "/"  # path to save the figure
    fig.savefig(path + figure_name, dpi=fig.dpi)
