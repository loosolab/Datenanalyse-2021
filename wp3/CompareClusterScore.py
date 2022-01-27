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
    parser.add_argument('Files', nargs="*", help='input files')
    parser.add_argument('-n', '--normalize', dest='Norm', default='Sum', help='Method of normalization of the data.')
    parser.add_argument('-o', '--outputName', dest='Custom_filename', default="ClusterComparison", help='Sets a custom name for saving the output files')
    parser.add_argument('-z', '--zscore', dest='Z', default=True, action='store_false', help='If used, the output will contain the absolute values (or those normlaized) and no z-Scores')
    args = parser.parse_args()
    
    return args.Files, args.Norm, args.Custom_filename, args.Z

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
    
    for i, element in enumerate(fileNames):
        df[element] =  df[element]*(100/normalizationFactors[i])
    
   # print(df)
    
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
    
    #print(dfZ)
    return dfZ

def main():
    possibleNormMethods = ["Sum","None"]
    
    inputFiles, normMethod, outputName, Z = cliParser()
    #outfile = open("ClusterComparison.tsv","a")
    
    #preparing a list of lists with columns for the TFs and each cluster with the scores for them
    scores = {"TF": []}
    transcriptionFactors = {}
    fileNames = []
    
    for element in inputFiles:
        file = open(element,"r")
        filename = os.path.basename(file.name)
        #print(filename)
        firstPartOfFilename, garbage = filename.split(".",1)
        columnName = firstPartOfFilename+"_Score"
        fileNames.append(columnName)
        scores[columnName] = []
        file.close()

    #print(scores)
    
    #getting the scores for each TF from each cluster, saving them in the list that was created beforehand
    for index,file in enumerate(inputFiles):

        currentFile = open(file, "r")
        columnName = fileNames[index]
        
        for currentIndex,line in enumerate(currentFile):
            #print(currentIndex,line)
            if line.startswith("output_prefix"):
                pass
            else:
                columns = line.split("\t")
                
                for element in columns:
                    try:
                        columns.remove("")
                    except:
                        pass
                
                
                if columns[1] in transcriptionFactors.keys():
                    transcriptionFactors[columns[1]][index] = columns[5]
                else:
                    transcriptionFactors[columns[1]] = []
                    for file in inputFiles:
                        transcriptionFactors[columns[1]].append(0)
                transcriptionFactors[columns[1]][index] = columns[5]
           
        currentFile.close()
                
    for key in transcriptionFactors.keys():
        scores["TF"].append(key)
        
        for index, cluster in enumerate(fileNames):
            scores[cluster].append(transcriptionFactors[key][index])

    
    #print(scores)    
    #print(transcriptionFactors)
    #checks whether normalization is wanted and which method is supposed to be used
    #normalizes the scores 
    #adjustedScores
    
    if normMethod not in possibleNormMethods:
        print("Normalization method not available, defaulting to 'Sum'")
        dfNormScores = normalizeToSum(scores, fileNames)
    else:    
        if normMethod == "Sum":
            dfNormScores = normalizeToSum(scores, fileNames)
        elif normMethod == "None":
            dfNormScores = pd.DataFrame(scores)
            dfNormScores.set_index("TF")

    if Z:        
        dfZScores = toZScore(dfNormScores, fileNames)  
    else:
        dfZScores = dfNormScores
            
    
    dfZScores.to_csv(outputName + '.tsv', sep="\t")
        

 
if __name__ == '__main__':
    main()
