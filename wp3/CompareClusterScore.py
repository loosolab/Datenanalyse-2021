# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 20:09:09 2021

@author: Moritz Hobein
"""

import argparse
import os.path
import pandas as pd
import numpy as np

#parser for using the tool via command-line interface
def cliParser():
    parser = argparse.ArgumentParser(description='Graphic Cluster Footprint Score Comparison')
    parser.add_argument('File', nargs=1, help='input file, bindetect_results.txt')
    parser.add_argument('-n', '--normalize', dest='Norm', default='Sum', help='Method of normalization of the data.')
    parser.add_argument('-o', '--outputName', dest='Custom_filename', default="ClusterComparison", help='Sets a custom name for saving the output files')
    parser.add_argument('-z', '--ZScore', dest='Z', default=True, action='store_false', help='if used, skips Z Scores calculation for each TF')
    parser.add_argument('-s', '--shortNames', dest='shortNames', default=False, action='store_true', help='If used, saves the short Name of each TF instead of the output_prefix')
    args = parser.parse_args()
    
    return args.File, args.Norm, args.Custom_filename, args.Z, args.shortNames

#normalizing the scores to the same sum
def normalizeToSum(scores, fileNames):
     
    normalizationFactors = []
    index = 0
    convertDict = {}
    
    for element in scores:
        if element != "TF":
            normalizationFactors.append(0)
            
            for value in scores[element]:
                normalizationFactors[index] = normalizationFactors[index] + float(value)
            index = index + 1  
    
    # for element in normalizationFactors:
    #     element = 1/element
    # print(normalizationFactors)
    
    df = pd.DataFrame(scores)
    df.set_index('TF', inplace=True)
    
    df = df.astype(float)
    #print(df)
    for i, element in enumerate(fileNames):
        df[element] =  df[element]*(100/normalizationFactors[i])
    
    #print(df)
    
    return df
       
#TODO
def toZScore(dfNorm, fileNames):
    deviations = []
    means = []
 
    dfZ_transposed = dfNorm.transpose()
    
    for column in dfZ_transposed:
        means.append(dfZ_transposed[column].mean())
        deviations.append(dfZ_transposed[column].std())
        mean = dfZ_transposed[column].mean()
        deviation = dfZ_transposed[column].std()
        dfZ_transposed[column] = (dfZ_transposed[column]-mean)/deviation
        
    #print(means)
    #print(deviations)
    #print(dfZ_transposed)
    
    dfZ = dfZ_transposed.transpose()
    
    print(dfZ)
    return dfZ
    

def main():
    possibleNormMethods = ["Sum","None"]
    
    inputFile, normMethod, outputName, Z, shortNames = cliParser()
    #outfile = open("ClusterComparison.tsv","a")
    
    #preparing a list of lists with columns for the TFs and each cluster with the scores for them
    scores = {"TF": []}
    transcriptionFactors = {}
    clusterNames = []
    clusterCol = []
    
    if shortNames:
        nameCol = 1
    else:
        nameCol = 0
    
    
    #inputData = pd.read_table(inputFile[0], sep="\t")
    file = open(inputFile[0], "r")
    
    for index, line in enumerate(file):
        if index == 0:
            noN, garbage = line.split("\n")
            #print(noN)
            elements = noN.split("\t")
            #print(elements)
            
            for i, col in enumerate(elements):
                if col.endswith("_mean_score"):
                    name = col.split("_mean_score")
                    clusterNames.append(name[0])
                    clusterCol.append(i)
                    scores[name[0]] = []
                    

    #print(scores)
    
    #getting the scores for each TF from each cluster, saving them in the list that was created beforehand
    #for index,file in enumerate(inputFiles):

        #currentFile = open(file, "r")
        #columnName = fileNames[index]
    file.seek(0)   
    for index ,line in enumerate(file):
        #print(currentIndex,line)
        if line.startswith("output_prefix"):
            pass
        else:
            noN = line.split("\n")
            columns = noN[0].split("\t")
            print(columns)   
                # for element in columns:
                #     try:
                #         columns.remove("")
                #     except:
                #         pass
            
            for i, value in enumerate(clusterCol):
                print(clusterCol)
                if columns[nameCol] in transcriptionFactors.keys():
                    transcriptionFactors[columns[nameCol]][i] = columns[value]
                else:
                    transcriptionFactors[columns[nameCol]] = []
                    for cluster in clusterNames:
                        transcriptionFactors[columns[nameCol]].append(0)
                    print(transcriptionFactors)
                    transcriptionFactors[columns[nameCol]][i] = columns[value]
               
    file.close()
                
    for key in transcriptionFactors.keys():
        scores["TF"].append(key)
        
        for index, cluster in enumerate(clusterNames):
            scores[cluster].append(transcriptionFactors[key][index])
    
    
    #print(clusterNames)
    
    print(scores)    
    #print(transcriptionFactors)
    #checks whether normalization is wanted and which method is supposed to be used
    #normalizes the scores 
    #adjustedScores
    
    if normMethod not in possibleNormMethods:
        print("Normalization method not available, defaulting to 'Sum'")
        dfNormScores = normalizeToSum(scores, clusterNames)
    else:    
        if normMethod == "Sum":
            dfNormScores = normalizeToSum(scores, clusterNames)
        elif normMethod == "None":
            dfNormScores = pd.DataFrame(scores)
            dfNormScores.set_index("TF")

    if Z:        
        dfZScores = toZScore(dfNormScores, clusterNames)  
    else:
        dfZScores = dfNormScores
            
    
    dfZScores.to_csv(outputName + '.tsv', sep="\t")
        

 
if __name__ == '__main__':
    main()