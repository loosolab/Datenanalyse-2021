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
    parser.add_argument('-f', '--families', dest='families', default=False, action='store_true', help='If used, output will group TFs by families')
    args = parser.parse_args()
    
    return args.families

def main():
    
    Families = cliParser()
    inputFile = "ClusterComparison.tsv"
    #inputFile = "comp.tsv"

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
    
    
    #print(clusters)
    
    collapseCols = {}
    
    dfT = df.transpose()
    
    if Families: 
        with open("TF_families.txt", "r") as file:
            for col in dfT:
                #print(dfT)
                for line in file:
                    family, tfactors = line.split("\t")
                    tfactorList = tfactors.split(",")
                    if col in tfactorList:
                            print(1)
                            collapseCols[col] = family
         
        dfFamilies = dfT
        print(dfFamilies)
        dfFamilies = dfFamilies.groupby(collapseCols, axis = 1).mean()
        print(dfFamilies)
        dfFamilies = dfFamilies.transpose()
        print(dfFamilies)
        
        
        df = dfFamilies
            
    print(collapseCols)    
    print(df)   
    #Quantilebestimmung und Extraktion der relevanten TF 
    quants = df.quantile(.95)
    
    
    for col in df:
        clusters[col] = []
    #print(quants["cluster9_Score"])
    
    for col in df:
        for index in df.index:
           if df.at[index, col] >= quants[col]:
               clusters[col].append(index)
    
    #print(clusters)
    important_TFs = pd.DataFrame(clusters)
    print(important_TFs)
    
    important_TFs.to_csv("importantTFsPerCluster.tsv" + '.tsv', sep="\t")

if __name__ == '__main__':
    main()