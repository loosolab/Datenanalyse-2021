#!/usr/bin/env python3

import os
import re
import pandas as pd
import plotly.express as px

print("###GENE SET ANALYSIS###")
print("Start gene set analysis")

## help method
# eliminates the \n at the end of the gene name
def extract_gene(gene_string: str):
    gene = gene_string[:-1]
    return gene


## analysis methods
# checks how many genes 2k has more than 1k
def compare_1k_2k(motif_1k: str, motif_2k: str, tissue_motif: str, cell_type_motif: str):
    # get motif name, tissue name and cell type name
    motif_name = ''
    tissue_name = tissue_motif
    cell_type_name = cell_type_motif
    # initialize counters
    counter_1k = 0
    counter_2k = 0
    ## read files
    # read file 1k
    with open(motif_1k) as file:
        for line in file:
            # check if line has a '#' in it
            if "#" in line:
                motif_name = line[1:-4]
                continue
            counter_1k += 1
    # read file 2k
    with open(motif_2k) as file:
        for line in file:
            # check if line has a '#' in it
            if "#" in line:
                continue
            counter_2k += 1
    # calculate the difference between the counters
    difference_2k_1k = counter_2k - counter_1k
    # calculate the difference in percent [%]
    difference_2k_1k_percent = difference_2k_1k / counter_2k * 100
    # returns the tissue name, the cell type name, the motif name, the counters and the difference as list
    return [tissue_name, cell_type_name, motif_name, counter_2k, counter_1k, difference_2k_1k, difference_2k_1k_percent]

# tests whether certain genes occur more than once
def identify_special_genes(motif_1k_all: str, motif_2k_all: str, tissue_motif: str, cell_type_motif: str):
    # get motif name, tissue name and cell type name
    motif_name = ''
    tissue_name = tissue_motif
    cell_type_name = cell_type_motif
    # initialize counters
    counter_1k_all = 0
    counter_2k_all = 0
    # initialize dictionaries
    dict_1k_all = {}
    dict_2k_all = {}
    # initialize ranking lists
    ranking_list_1k = []
    ranking_list_2k = []
    ## read files
    # read file 1k all
    with open(motif_1k_all) as file:
        for line in file:
            # check if line has a '#' in it
            if "#" in line:
                motif_name = line[1:-4]
                continue
            counter_1k_all += 1
            # extract gene name
            gene = extract_gene(line)
            # check if gene is already in dictionary
            # if not: add gene to dictionary
            if gene not in dict_1k_all.keys():
                dict_1k_all[gene] = 1
            # if in dict: increase counter
            else:
                val = dict_1k_all[gene] + 1
                dict_1k_all[gene] = val
    # read file 2k all
    with open(motif_2k_all) as file:
        for line in file:
            # check if line has a '#' in it
            if "#" in line:
                continue
            counter_2k_all += 1
            # extract gene name
            gene = extract_gene(line)
            # check if gene is already in dictionary
            # if not: add gene to dictionary
            if gene not in dict_2k_all.keys():
                dict_2k_all[gene] = 1
            # if in dict: increase counter
            else:
                val = dict_2k_all[gene] + 1
                dict_2k_all[gene] = val
    # calculate the difference between the counters
    difference_2k_all_1k_all = counter_2k_all - counter_1k_all
    # generate ranking tables
    dict_1k_all = (sorted(dict_1k_all.items(), key=lambda item: item[1], reverse=True))
    dict_2k_all = (sorted(dict_2k_all.items(), key=lambda item: item[1], reverse=True))
    # generate ranking lists by output_genes
    for ranked_gene in dict_1k_all:
        if ranked_gene[1] > 1:
            ranking_list_1k.append(ranked_gene)
    
    for ranked_gene in dict_2k_all:
        if ranked_gene[1] > 1:
            ranking_list_2k.append(ranked_gene)
            
    # returns the tissue name, the cell type name, the motif name, the counters and the difference as list and the ranked gene lists (2k + 1k) of the genes that occur more than once
    return [tissue_name, cell_type_name , motif_name, counter_2k_all, counter_1k_all, difference_2k_all_1k_all], ranking_list_2k, ranking_list_1k
    
# tests whether a motif has a similarity to a known transcription factor or other motifs
def correlation_to_TF_or_motif(motif: str , TF_or_motif: str, tissue_motif: str, cell_type_motif: str, tissue_TF_or_motif: str, cell_type_TF_or_motif: str): 
    # initialize motif name
    motif_name = ''
    # initialize naming list for TFs or motifs -> depends on analysis
    TF_or_motif_name_list = []
    # get tissue and cell type names for the analyzed files
    tissue_motif_name = tissue_motif
    cell_type_motif_name = cell_type_motif
    tissue_TF_or_motif_name = tissue_TF_or_motif
    cell_type_TF_or_motif_name = cell_type_TF_or_motif
    # initialize counters: for motif and for TFs or motifs
    counter_motif = 0
    counter_TF_or_motif_list = []
    # initialize lists for counters of matching genes, genes of the motif, genes of the TFs or motifs and the similarity
    counter_matching_genes_list = []
    motif_gene_list = []
    TFs_or_motifs_gene_list = []
    correlated_TFs_or_motifs_similarity_list = []
    ## read files
    # read motif
    with open(motif) as file:
        for line in file:
            # check if line has a '#' in it
            if "#" in line:
                # extract motif name
                motif_name = line[1:-4]
                continue
            counter_motif += 1
            # extract gene names
            gene = extract_gene(line)
            motif_gene_list.append(gene)
    # read TF_or_motif to count lines
    # initialize counters for lines in file and current line
    line_counter = 0
    current_line = 0
    with open(TF_or_motif) as file:
        for line in file:
            line_counter += 1
    # read TF_or_motif again
    with open(TF_or_motif) as file:
        # initialize help list and counter
        TF_or_motif_gene_list = []
        counter_TF_or_motif = 0
        for line in file:
            current_line += 1
            # check if line has a '#' in it
            if "#" in line:
                # extract TF or motif name and add to list
                TF_or_motif_name = line[1:-4]
                TF_or_motif_name_list.append(TF_or_motif_name)                
            # check if line is empty
            # if not empty: increase counter, extract gene name, append gene name to list
            elif(line != '\n'):
                counter_TF_or_motif += 1
                gene = extract_gene(line)
                TF_or_motif_gene_list.append(gene)
                # check if line is last line
                # if last line: append help gene list to gene list, add help counter to counter list
                if(current_line == line_counter):
                    TFs_or_motifs_gene_list.append(TF_or_motif_gene_list)
                    counter_TF_or_motif_list.append(counter_TF_or_motif)
            # if empty: add help gene list to gene list, add help counter to counter list
            # reset help list and counter
            elif(line == '\n'):
                TFs_or_motifs_gene_list.append(TF_or_motif_gene_list)
                counter_TF_or_motif_list.append(counter_TF_or_motif)
                TF_or_motif_gene_list = []
                counter_TF_or_motif = 0
                
    ## calculations
    # calculate the matching genes -> fill list: counter_matching_genes_list
    for TF_or_motif in TFs_or_motifs_gene_list:
        counter_matching_genes = 0
        for gene_TF_or_motif, motif_gene in zip(TF_or_motif, motif_gene_list):
            if(gene_TF_or_motif == motif_gene):
                counter_matching_genes += 1
        counter_matching_genes_list.append(counter_matching_genes)
                    
    # calculate the correlated TFs for the motif
    # -> fill list: correlated_TFs_similarity_list
    index_counter = 0
    # check matching counters
    for matching_genes in counter_matching_genes_list:
        # matching counter not 0: calculate similarity
        if(matching_genes != 0):
            similarity = (matching_genes / counter_motif) * 100
            correlated_TFs_or_motifs_similarity_list.append(similarity)
            index_counter += 1
        # matching counter is 0: similarity = 0
        else:
            similarity = 0
            correlated_TFs_or_motifs_similarity_list.append(similarity)
    # returns tissue name of the motif, cell type name of the motif, motif name, counted genes of motif, tissue name of the other motif or TF, cell type name of the other motif or TF
    # also returns some lists: naming list of the motifs or TFs, list with counted genes of motifs or TFs, list with the counter of matching genes to the motif, list with the similarity value to the motif [%]
    return [tissue_motif_name, cell_type_motif_name, motif_name, counter_motif, tissue_TF_or_motif_name, cell_type_TF_or_motif_name], TF_or_motif_name_list, counter_TF_or_motif_list, counter_matching_genes_list, correlated_TFs_or_motifs_similarity_list 


## get file paths and store them in lists
# initialize lists for file paths
# stores the paths of the motifs to be compared -> one list in list = one motif pair (motifx_1k + motifx_2k)
# compares the same motif itself with maybe different gene sets
file_paths_compare_1k_2k = []
# stores the paths of the motifs to be compared -> one list in list = one motif pair (motifx_1k_all + motifx_2k_all)
# compares the same motif itself with a look at the multiple occurring genes
file_paths_identify_special_genes = []
# stores the paths of the motifs and TFs to be compared -> one list in list = one motif TF pair (motifx_2k + TFy_2k)
# compares the motif with a set of TFs (TF set = per tissue) to identify correlations
file_paths_correlation_to_TF = []
# stores the paths of the motifs to be compared -> one list in list = one motif pair (motifx_2k + motify_2k)
# compares the motif with a set of motifs (motif set = per cell type) to identify correlations
file_paths_correlation_to_motif = []

# get current working directory
path = os.getcwd()

# change path to runs
os.chdir('../runs')
runs_path = os.getcwd()

# create help lists
# help list for file_paths_compare_1k_2k
motif_list_compare_1k_2k = []
# help list for file_paths_identify_special_genes
motif_list_identify_special_genes = []
# help lists for file_paths_correlation_to_TF and file_paths_correlation_to_motif
# store all TF paths
TFs = []
# store all cell_type motif paths
motifs_cell_type = []
# store all motifs
motifs_for_TFs_and_motifs_cell_type = []

# collect paths to files   
print("Collect file paths")
for tissue in os.listdir(runs_path):
    t = os.path.join(runs_path, tissue)
    if(os.path.isdir(t)):
        #print(tissue)
        os.chdir(t)
        tissue_path = os.getcwd()
        #print(tissue_path)
        
        # TF gene sets are stored here -> files
        for filename_TF in os.listdir(tissue_path):
            TF = os.path.join(tissue_path, filename_TF)
            # checking if it is a file
            if(os.path.isfile(TF)):
                #print(TF)
                if re.match(".*[2]k.txt", filename_TF):
                    TFs.append(TF)
                
        for cell_type in os.listdir(tissue_path):
            ct = os.path.join(tissue_path, cell_type)
            if(os.path.isdir(ct)):
                #print(cell_type)
                os.chdir(ct)
                cell_type_path = os.getcwd()
                #print(cell_type_path)
                
                # cell type gene sets are stored here -> files
                os.chdir('./annotation')
                annotation_path = os.getcwd()
                #print(annotation_path)
                for filename_cell_type in os.listdir(annotation_path):
                    # ctm = cell type motifs
                    ctm = os.path.join(annotation_path, filename_cell_type)
                    # checking if it is a file
                    if(os.path.isfile(ctm)):
                        #print(ctm)
                        if re.match(".*[2]k.txt", filename_cell_type):
                            motifs_cell_type.append(ctm)
                        
                # motif gene sets are stored here -> files
                os.chdir('./gene_sets_motifs')
                gene_sets_motifs_path = os.getcwd()
                #print(gene_sets_motifs_path)
                # initialize help list for compare_1k_2k
                motif_list = []
                # initialize help list for identify_special_genes
                motif_list_all = []
                for filename_motif in os.listdir(gene_sets_motifs_path):
                    # m = cell type motifs
                    m = os.path.join(gene_sets_motifs_path, filename_motif)
                    # checking if it is a file
                    if(os.path.isfile(m)):
                        #print(m)
                        if re.match(".*[1-2]k.txt", filename_motif):
                            motif_list.append(m)
                            if re.match(".*[2]k.txt", filename_motif):
                                motifs_for_TFs_and_motifs_cell_type.append(m)
                        else:
                            motif_list_all.append(m)
                # add help lists to other help lists            
                motif_list_compare_1k_2k.append(motif_list)
                motif_list_identify_special_genes.append(motif_list_all)
       
# fill path lists 
# fill file_paths_compare_1k_2k
for motifs in motif_list_compare_1k_2k:
    for i in range(0, len(motifs), 2):
        compare_list = [motifs[i], motifs[i+1]]
        file_paths_compare_1k_2k.append(compare_list)
        
# fill file_paths_identify_special_genes
for motifs in motif_list_identify_special_genes:
    for i in range(0, len(motifs), 2):
        compare_list = [motifs[i], motifs[i+1]]
        file_paths_identify_special_genes.append(compare_list)
        
# fill file_paths_correlation_to_TF
for motif in motifs_for_TFs_and_motifs_cell_type:
    compare_list = []
    for TF in TFs:
        # identify tissue
        # motif (tissue)
        string_motif = motif.split('/')
        tissue_motif = string_motif[-5]
        # TF (tissue) 
        string_TF = TF.split('/')
        tissue_TF = string_TF[-2]
        # test whether motif and transcription factor belong to the same tissue
        if(tissue_motif == tissue_TF):
            compare_list = [motif, TF]
            file_paths_correlation_to_TF.append(compare_list)

# fill file_paths_correlation_to_motif
for motif in motifs_for_TFs_and_motifs_cell_type:
    compare_list = []
    for cell_type in motifs_cell_type:
        compare_list = [motif, cell_type]
        file_paths_correlation_to_motif.append(compare_list)

# delete help lists
del motif_list_compare_1k_2k
del motif_list_identify_special_genes
del TFs
del motifs_cell_type
del motifs_for_TFs_and_motifs_cell_type
                        
# change to origin directory
os.chdir(path)


## calculations
print("Start calculations")
# output lists
compare_1k_2k_results = []
identify_special_genes_results = []
correlation_to_TF_results = []
correlation_to_motif_results = []

# run motif analysis for differences in gene sets spaced 1000 and 2000 base pairs apart, respectively
print("1. Compare difference in gene sets")
for elem in file_paths_compare_1k_2k:
    # extract tissue and cell type name
    tissue_cell_type = elem[0].split('/')
    tissue = tissue_cell_type[-5]
    cell_type = tissue_cell_type[-4]
    result = compare_1k_2k(elem[0], elem[1], tissue, cell_type)
    compare_1k_2k_results.append(result)

# run motif analysis to identify genes that have been annotated multiple times
print("2. Identify special genes")
for elem in file_paths_identify_special_genes:
    # extract tissue and cell type name
    tissue_cell_type = elem[0].split('/')
    tissue = tissue_cell_type[-5]
    cell_type = tissue_cell_type[-4]
    result = identify_special_genes(elem[0], elem[1], tissue, cell_type)
    # create lists to eliminate the nested list
    list_1 = result[0]
    list_2 = result[1]
    list_3 = result[2]
    # unpack the tuples in list_2 and list_3
    for tup_2, tup_3 in zip(list_2, list_3):
        # create list to save the correctly modified output
        ready_list = []
        for element in list_1:
            ready_list.append(element)
        ready_list.append(tup_2[0])
        ready_list.append(tup_2[1])
        ready_list.append(tup_3[0])
        ready_list.append(tup_3[1])
        identify_special_genes_results.append(ready_list)
    
# run motif analysis to identify similar TFs
print("3. Identify similar TFs")
for elem in file_paths_correlation_to_TF:
    # extract tissue and cell type name for motif and TF
    tissue_cell_type = elem[0].split('/')
    tissue_motif = tissue_cell_type[-5]
    cell_type_motif = tissue_cell_type[-4]
    tissue_cell_type_TF = elem[1].split('/')
    tissue_TF = tissue_cell_type_TF[-2]
    # TF does not have an associated cell type -> 'NA'
    cell_type_TF = 'NA'
    result = correlation_to_TF_or_motif(elem[0], elem[1], tissue_motif, cell_type_motif, tissue_TF, cell_type_TF)
    # create lists to eliminate the nested list
    list_1 = result[0]
    list_2 = result[1]
    list_3 = result[2]
    list_4 = result[3]
    list_5 = result[4]
    for i in range(len(list_2)):
        # create list to save the correctly modified output
        ready_list = []
        for element in list_1:
            ready_list.append(element)
        ready_list.append(list_2[i])
        ready_list.append(list_3[i])
        ready_list.append(list_4[i])
        ready_list.append(list_5[i])
        correlation_to_TF_results.append(ready_list)

# run motif analysis to identify similar and same motifs
print("4. Identify similar and same motifs")
for elem in file_paths_correlation_to_motif:
    # extract tissue and cell type name for motif and other motifs
    tissue_cell_type = elem[0].split('/')
    tissue_motif = tissue_cell_type[-5]
    cell_type_motif = tissue_cell_type[-4]
    tissue_cell_type_motifs = elem[1].split('/')
    tissue_motifs = tissue_cell_type_motifs[-4]
    cell_type_motifs = tissue_cell_type_motifs[-3]
    result = correlation_to_TF_or_motif(elem[0], elem[1], tissue_motif, cell_type_motif, tissue_motifs, cell_type_motifs)
    # create lists to eliminate the nested list
    list_1 = result[0]
    list_2 = result[1]
    list_3 = result[2]
    list_4 = result[3]
    list_5 = result[4]
    for i in range(len(list_2)):
        # create list to save the correctly modified output
        ready_list = []
        for element in list_1:
            ready_list.append(element)
        ready_list.append(list_2[i])
        ready_list.append(list_3[i])
        ready_list.append(list_4[i])
        ready_list.append(list_5[i])
        correlation_to_motif_results.append(ready_list)
        

## prepare Dataframes
df_compare_1k_2k = pd.DataFrame(compare_1k_2k_results, columns=['tissue','cell_type','motif','number_of_genes_2k','number_of_genes_1k','difference_number_of_genes', 'difference_in_percent[%]'])
df_identify_special_genes = pd.DataFrame(identify_special_genes_results, columns=['tissue', 'cell_type', 'motif', 'number_of_genes_2k_all', 'number_of_genes_1k_all', 'difference_number_of_genes', 'gene_2k', 'counter_gene_2k', 'gene_1k', 'counter_gene_1k'])
df_correlation_to_TF = pd.DataFrame(correlation_to_TF_results, columns=['tissue', 'cell_type', 'motif', 'number_of_genes_2k', 'tissue_TF', 'cell_type_TF', 'TF', 'number_of_genes_TF_2k', 'number_of_matching_genes', 'similarity[%]'])
df_correlation_to_motif = pd.DataFrame(correlation_to_motif_results, columns=['tissue', 'cell_type', 'motif', 'number_of_genes_2k', 'tissue_other_motif', 'cell_type_other_motif', 'other_motif', 'number_of_genes_other_motif_2k', 'number_of_matching_genes', 'similarity[%]'])

# delete output lists
del compare_1k_2k_results
del identify_special_genes_results
del correlation_to_TF_results
del correlation_to_motif_results


## write dataframes to csv
# create output folder if not existing
outdir = '../gene_set_analysis'
if not os.path.exists(outdir):
    os.mkdir(outdir)
# creating outputs
# output compare_1k_2k
outname_compare_1k_2k = 'compare_1k_2k.csv'
fullname_compare_1k_2k = os.path.join(outdir, outname_compare_1k_2k)
df_compare_1k_2k.to_csv(fullname_compare_1k_2k)

# output identify_special_genes
outname_identify_special_genes = 'identify_special_genes.csv'
fullname_identify_special_genes = os.path.join(outdir, outname_identify_special_genes)
df_identify_special_genes.to_csv(fullname_identify_special_genes)
# reduce with cutoff
cutoff_counter_genes = 15
df_identify_special_genes_red = df_identify_special_genes[df_identify_special_genes['difference_number_of_genes'] > cutoff_counter_genes].copy()
# write to csv
outname_identify_special_genes_red = 'identify_special_genes_red.csv'
fullname_identify_special_genes_red = os.path.join(outdir, outname_identify_special_genes_red)
df_identify_special_genes_red.to_csv(fullname_identify_special_genes_red)

# output correlation to TF
outname_correlation_to_TF = 'correlation_to_TF.csv'
fullname_correlation_to_TF = os.path.join(outdir, outname_correlation_to_TF)
df_correlation_to_TF.to_csv(fullname_correlation_to_TF)
# reduce with cutoff
cutoff_TF = 1
df_correlation_to_TF_red = df_correlation_to_TF[df_correlation_to_TF['similarity[%]'] > cutoff_TF].copy()
# write to csv
outname_correlation_to_TF_red = 'correlation_to_TF_red.csv'
fullname_correlation_to_TF_red = os.path.join(outdir, outname_correlation_to_TF_red)
df_correlation_to_TF_red.to_csv(fullname_correlation_to_TF_red)

# output correlation to motif
outname_correlation_to_motif = 'correlation_to_motif.csv'
fullname_correlation_to_motif = os.path.join(outdir, outname_correlation_to_motif)
df_correlation_to_motif.to_csv(fullname_correlation_to_motif)
# reduce with cutoff_1_motif and cutoff_2_motif
cutoff_1_motif = 100
cutoff_2_motif = 1
df_correlation_to_motif_red_1 = df_correlation_to_motif[df_correlation_to_motif['similarity[%]'] < cutoff_1_motif].copy()
df_correlation_to_motif_red_2 = df_correlation_to_motif_red_1[df_correlation_to_motif_red_1['similarity[%]'] > cutoff_2_motif].copy()
# write to csv
outname_correlation_to_motif_red = 'correlation_to_motif_red.csv'
fullname_correlation_to_motif_red = os.path.join(outdir, outname_correlation_to_motif_red)
df_correlation_to_motif_red_2.to_csv(fullname_correlation_to_motif_red)


### plotting
colorlist = ["aliceblue", "antiquewhite", "aqua", "aquamarine", "azure", "beige", "bisque", "black", "blanchedalmond", "blue", "blueviolet", "brown", "burlywood", "cadetblue",
             "chartreuse", "chocolate", "coral", "cornflowerblue", "cornsilk", "crimson", "cyan", "darkblue", "darkcyan", "darkgoldenrod", "darkgray", "darkgrey", "darkgreen",
             "darkkhaki", "darkmagenta", "darkolivegreen", "darkorange", "darkorchid", "darkred", "darksalmon", "darkseagreen", "darkslateblue", "darkslategray", "darkslategrey",
             "darkturquoise", "darkviolet", "deeppink", "deepskyblue", "dimgray", "dimgrey", "dodgerblue", "firebrick", "floralwhite", "forestgreen", "fuchsia", "gainsboro",
             "ghostwhite", "gold", "goldenrod", "gray", "grey", "green", "greenyellow", "honeydew", "hotpink", "indianred", "indigo", "ivory", "khaki", "lavender", "lavenderblush", "lawngreen",
             "lemonchiffon", "lightblue", "lightcoral", "lightcyan", "lightgoldenrodyellow", "lightgray", "lightgrey", "lightgreen", "lightpink", "lightsalmon", "lightseagreen",
             "lightskyblue", "lightslategray", "lightslategrey", "lightsteelblue", "lightyellow", "lime", "limegreen", "linen", "magenta", "maroon", "mediumaquamarine",
             "mediumblue", "mediumorchid", "mediumpurple", "mediumseagreen", "mediumslateblue", "mediumspringgreen", "mediumturquoise", "mediumvioletred", "midnightblue",
             "mintcream", "mistyrose", "moccasin", "navajowhite", "navy", "oldlace", "olive", "olivedrab", "orange", "orangered", "orchid", "palegoldenrod", "palegreen", "paleturquoise",
             "palevioletred", "papayawhip", "peachpuff", "peru", "pink", "plum", "powderblue", "purple", "red", "rosybrown", "royalblue", "saddlebrown", "salmon", "sandybrown",
             "seagreen", "seashell", "sienna", "silver", "skyblue", "slateblue", "slategray", "slategrey", "snow", "springgreen", "steelblue", "tan", "teal", "thistle", "tomato", "turquoise",
             "violet", "wheat", "white", "whitesmoke", "yellow", "yellowgreen"]

## create percentage bar plots for compare_1k_2k 
# create unique Motif and Cell_Type name and add to Data Frame
df_compare_1k_2k['Motif'] = df_compare_1k_2k['tissue'] + "_" + df_compare_1k_2k['cell_type'] + "_" + df_compare_1k_2k['motif']
df_compare_1k_2k['Cell_Type'] = df_compare_1k_2k['tissue'] + "_" + df_compare_1k_2k['cell_type']

## Bar plots
# Tissue
fig_1 = px.bar(df_compare_1k_2k, x="Motif", y="difference_in_percent[%]", color="tissue", title="Percentage Deviation of Genes from 2k to 1k per Tissue", color_discrete_sequence=colorlist)
fig_1.update_layout(barmode='stack', xaxis={'categoryorder':'total descending'})
fig_1.update_traces(marker_line_width=0)
# save
fig_1.write_html(f"../gene_set_analysis/compare_1k_2k_bar_Tis.html")
# Cell Type
fig_1 = px.bar(df_compare_1k_2k, x="Motif", y="difference_in_percent[%]", color="Cell_Type", title="Percentage Deviation of Genes from 2k to 1k per Cell Type", color_discrete_sequence=colorlist)
fig_1.update_layout(barmode='stack', xaxis={'categoryorder':'total descending'})
fig_1.update_traces(marker_line_width=0)
# save
fig_1.write_html(f"../gene_set_analysis/compare_1k_2k_bar_CT.html")

## create special gene count bar plots for identify_special_genes
# cutoff = 15
# create unique Motif and Cell_Type name and add to Data Frame
df_identify_special_genes_red['Motif'] = df_identify_special_genes_red['tissue'] + "_" + df_identify_special_genes_red['cell_type'] + "_" + df_identify_special_genes_red['motif']
df_identify_special_genes_red['Cell_Type'] = df_identify_special_genes_red['tissue'] + "_" + df_identify_special_genes_red['cell_type']

## Bar plots
# All motifs
fig_2 = px.bar(df_identify_special_genes_red, x="gene_2k", y="counter_gene_2k", color="Motif", title="Counts of multiple annotated Genes across all Motifs at a cutoff of 15 for all Motifs", color_discrete_sequence=colorlist)
fig_2.update_layout(barmode='stack', xaxis={'categoryorder':'total descending'})
fig_2.update_traces(marker_line_width=0)
# save
fig_2.write_html(f"../gene_set_analysis/identify_special_genes_bar_Mot.html")
# Tissue
fig_2 = px.bar(df_identify_special_genes_red, x="gene_2k", y="counter_gene_2k", color="tissue", title="Counts of multiple annotated Genes across all Motifs at a cutoff of 15 for all Tissues", color_discrete_sequence=colorlist)
fig_2.update_layout(barmode='stack', xaxis={'categoryorder':'total descending'})
fig_2.update_traces(marker_line_width=0)
# save
fig_2.write_html(f"../gene_set_analysis/identify_special_genes_bar_Tis.html")
# Cell Type
fig_2 = px.bar(df_identify_special_genes_red, x="gene_2k", y="counter_gene_2k", color="Cell_Type", title="Counts of multiple annotated Genes across all Motifs at a cutoff of 15 for all Cell Types", color_discrete_sequence=colorlist)
fig_2.update_layout(barmode='stack', xaxis={'categoryorder':'total descending'})
fig_2.update_traces(marker_line_width=0)
# save
fig_2.write_html(f"../gene_set_analysis/identify_special_genes_bar_CT.html")

## create bar plot for correlation_to_TF
# cutoff = 1 %
# create motif names
df_correlation_to_TF_red['Motif'] = df_correlation_to_TF_red['tissue'] + "_" + df_correlation_to_TF_red['cell_type'] + "_" + df_correlation_to_TF_red['motif']

## Bar plot
# All motifs
fig_3 = px.bar(df_correlation_to_TF_red, x="Motif", y="similarity[%]", color="TF", title="Similarity between motifs and TFs", color_discrete_sequence=colorlist)
fig_3.update_layout(barmode='stack', xaxis={'categoryorder':'total descending'})
fig_3.update_traces(marker_line_width=0)
# save
fig_3.write_html(f"../gene_set_analysis/correlation_to_TF_bar_Mot.html")

## create bar plot for correlation_to_motif
# cutoff = 1 % and 100 % -> exclude the motifs that have matched themselves
# create motif names for the compared motifs
df_correlation_to_motif_red_2['Motif_1'] = df_correlation_to_motif_red_2['tissue'] + "_" + df_correlation_to_motif_red_2['cell_type'] + "_" + df_correlation_to_motif_red_2['motif']
df_correlation_to_motif_red_2['Motif_2'] = df_correlation_to_motif_red_2['tissue_other_motif'] + "_" + df_correlation_to_motif_red_2['cell_type_other_motif'] + "_" + df_correlation_to_motif_red_2['other_motif']

## Bar plot
# All motifs
fig_4 = px.bar(df_correlation_to_motif_red_2, x="Motif_2", y="similarity[%]", color="Motif_1", title="Similarity between motifs", color_discrete_sequence=colorlist)
fig_4.update_layout(barmode='stack', xaxis={'categoryorder':'total descending'})
fig_4.update_traces(marker_line_width=0)
# save
fig_4.write_html(f"../gene_set_analysis/correlation_to_motif_bar_Mot.html")


print("Done")
