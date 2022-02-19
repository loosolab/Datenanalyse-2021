#!/usr/bin/env python3
import matplotlib.pyplot as plt
import pandas as pd

# Initialize data lists for tissue, cell_type und motif_counts
data = pd.read_csv('/mnt/workspace_stud/allstud/wp5/runs/motif_counts.txt', sep='\t')

df = pd.DataFrame(data)

tissues = list(df.iloc[:, 0])
cell_types = list(df.iloc[:, 1])
motif_counts = list(df.iloc[:, 2])

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
    cell_type_x = cell_types[0:tissue_count]
    motif_counts_x = motif_counts[0:tissue_count]

    # remove used values from data lists
    tissues = tissues[tissue_count:]
    cell_types = cell_types[tissue_count:]
    motif_counts = motif_counts[tissue_count:]

    # change the figure size
    fig = plt.figure(figsize=(18, 14), dpi=300)

    # plot the data as bar plot
    plt.bar(cell_type_x, motif_counts_x, color='b')
    plt.xticks(rotation=90)
    plt.title(new_tissue)
    plt.xlabel("Cell Types", labelpad=15)
    plt.ylabel("Motif Counts")

    # show bar plot
    plt.show()
    
    # save plot
    figure_name = "motif_counts_" + new_tissue + "_barplot.png"
    path = "/mnt/workspace_stud/allstud/wp5/runs/" + new_tissue + "/"  # path to save the figure
    fig.savefig(path + figure_name, dpi=fig.dpi)
