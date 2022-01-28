# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 19:51:13 2022

@author: Moritz Hobein
"""
import pandas as pd
import argparse

#parser for using the tool via command-line interface
def cliParser():
    parser = argparse.ArgumentParser(description='Graphic Cluster Footprint Score Comparison')
    parser.add_argument('File', nargs=1, help='input file, tsv of binding scores')
    parser.add_argument('-n', '--normalize', dest='Norm', default='Sum', help='Method of normalization of the data.')
    parser.add_argument('-o', '--outputName', dest='Custom_filename', default="ClusterComparison", help='Sets a custom name for saving the output files')
    parser.add_argument('-z', '--zscore', dest='Z', default=True, action='store_false', help='If used, the output will contain the absolute values (or those normlaized) and no z-Scores')
    args = parser.parse_args()
    
    return args.File

def main():
    
    inputFile = "comp.tsv"

    with open(inputFile, "r") as file:
        
        headerCount = 0
        
        for line in file:
        
            if line.startswith("#") or line.startswith("!"):
               pass
            elif line.startswith("TF"):
               break
            headerCount += 1
               
    
    df = pd.read_csv(inputFile, sep="\t", header = headerCount)
    df.set_index('TF', inplace=True)

    clusters = {}
    
    for col in df:
        clusters[col] = []
    
    print(clusters)
    
    
    #Quantilebestimmung und Extraktion der relevanten TF 
    quants = df.quantile(.95)
    
    print(quants["cluster9_Score"])
    
    for col in df:
        for index in df.index:
           if df.at[index, col] >= quants[col]:
               clusters[col].append(index)
    
    print(clusters)
    important_TFs = pd.DataFrame(clusters)
    print(important_TFs)
    
    important_TFs.to_csv("importantTFsPerCluster.tsv" + '.tsv', sep="\t")

if __name__ == '__main__':
    main()