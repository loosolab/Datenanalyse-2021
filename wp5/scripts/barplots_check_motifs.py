#!/usr/bin/env python3
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Initialize data lists for tissue, cell_type und motifs
data = pd.read_csv('/mnt/workspace_stud/allstud/wp5/runs/check_new_motifs_all.tx', sep='\t')
tissue_names = pd.unique(data['tissue'])

# create plot for all tissues

sns.set_style("ticks")
# plot the data as bar plot
fig = sns.countplot(data=data, x="sequence", order=data["sequence"].value_counts().index)
plt.xticks(rotation=90)
plt.title("Motif Counts for all new Motifs (all tissues)")
plt.xlabel("New Motifs", labelpad=15)
sns.despine()

# show bar plot
plt.show()

# save plot
path = "/mnt/workspace_stud/allstud/wp5/runs/"  # path to save the figure
figure_name = "new_motif_counts_for_all_tissues_barplot.png"
fig.savefig(path + figure_name, dpi=fig.dpi)
plt.close()

# create plot per tissue
for tissue in tissue_names:
    # plot the data as bar plot
    tis_data = data[data['tissue'] == tissue]

    ## plot barplot per motif
    fig = sns.countplot(data=tis_data, x="sequence", order=tis_data["sequence"].value_counts().index)
    plt.xticks(rotation=90)
    plt.title("Motif Counts for all new Motifs (" + tissue + ")")
    plt.xlabel("New Motifs", labelpad=15)

    # show bar plot
    plt.show()

    # save plot
    figure_name = "new_motif_counts_for_" + tissue + "_barplot.png"
    path = "/mnt/workspace_stud/allstud/wp5/runs/" + tissue + "/"  # path to save the figure
    fig.savefig(path + figure_name, dpi=fig.dpi)
    plt.close()

    ## plot bar plot per cluster
    fig = sns.countplot(data=tis_data, x="cell_type", order=tis_data["cell_type"].value_counts().index)
    plt.xticks(rotation=90)
    plt.title("Cell type Counts for all new Motifs (" + tissue + ")")
    plt.xlabel("Cell types", labelpad=15)

    # show bar plot
    plt.show()

    # save plot
    figure_name = figure_name = "motif_counts_" + new_tissue + "_barplot.png"
    path = "/mnt/workspace_stud/allstud/wp5/runs/" + tissue + "/"  # path to save the figure
    fig.savefig(path + figure_name, dpi=fig.dpi)
    plt.close()
