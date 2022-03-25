# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 19:51:13 2022

@author: Moritz Hobein
"""
import pandas as pd
import argparse

#parser for using the tool via command-line interface
def cliParser():
    parser = argparse.ArgumentParser(description='Script to extract the transcription factors that define a cluster (by z score of TOBIAS footprinting score')
    parser.add_argument('File', nargs=1, help='input file, tsv of binding scores')
    parser.add_argument('-o', '--outputName', dest='Custom_filename', default="importantTFsPerCluster.tsv", help='Sets a custom name for saving the output files')
    parser.add_argument('-q', '--quantile', dest='Quantile', default=.95, help='Sets quantile of what the defining TFs are. Default is .95 (give out the top 5% of TFs by z score)')
    args = parser.parse_args()
    
    return args.File, args.Quantile, args.Custom_filename

def getDefiningTF(inputFile, quantile, outname):
    
    #counting the header lines to ignore them in data extraction later
    with open(inputFile, "r") as file:
        
        headerCount = 0
        
        for line in file:
        
            if line.startswith("#") or line.startswith("!"):
               pass
            elif line.startswith("TF"):
               break
            headerCount += 1
               
    #reading the input file whole ignoring the perviously established header
    df = pd.read_csv(inputFile, sep="\t", header = headerCount)
    df.set_index('TF', inplace=True)

    clusters = {}
    
    #extracting the top TFs by z score according to the set quantile
    quants = df.quantile(quantile)
    
    for col in df:
        clusters[col] = []

    
    for col in df:
        for index in df.index:
           if df.at[index, col] >= quants[col]:
               clusters[col].append(index)
    
    #creating the output
    important_TFs = pd.DataFrame(clusters)
    
    important_TFs.to_csv(outname, sep="\t")
    
def main():
    
    File, quantile, outname = cliParser()
    inputFile = File[0]

    getDefiningTF(inputFile, quantile, outname)

if __name__ == '__main__':
    main()
