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
    parser.add_argument('-f', '--families', dest='families', nargs=1, default=[""], help='Family cluster file to group TF families together')
    args = parser.parse_args()
    
    return args.File, args.Norm, args.Custom_filename, args.Z, args.shortNames, args.families

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
    
    #print(dfZ)
    return dfZ
    
def groupFamilies(familyFile, df):
   
    with open(familyFile, "r") as file:
        
        colnames = []
        
        for col in df:
            colnames.append(col)
        
        dfT = df.transpose()
        
        families = {}

        for line in file:
            #print(line)
            familyName, familyMembers = line.split("\t")
            familyMembers= familyMembers.split("\n")
            familyMembers = familyMembers[0].split(",")
            
            for index, member in enumerate(familyMembers):
                familyMembers[index] = member.upper()
            
            families[familyName] = familyMembers
        

        familyScores = {}
        
        for familyName in families.keys():
            #print(familyName)
            tmpScores = []
            memberCounter = 0
            
            
            for col in dfT:
                #print(col)
                #print(families[familyName])
                if col.upper() in families[familyName]:
                    if tmpScores == []:
                        tmpScores = dfT[col].tolist()
                        memberCounter += 1
                    else:
                        for i, score in enumerate(dfT[col].tolist()):
                            tmpScores[i] += score
                        memberCounter += 1
                
                    if memberCounter == len(families[familyName]):
                        break
                #print(memberCounter)
           # print(tmpScores)
            
            meanScores = []
            
            for score in tmpScores:
                meanScores.append(score/memberCounter)
            
            familyScores[familyName] = meanScores
        
        
        # for element in familyScores:
        #     if len(familyScores[element]) == 0:
        #         print(element,": ",familyScores[element],": ", df[element].tolist())

        #TODO so soll das nicht, warum muss ich welche rauswerfen, weil sie leer sind????
        finale = {}
        
        for element in familyScores:
            if len(familyScores[element]) != 0:
                finale[element] = familyScores[element]
            else:
                finale[element] = []
                
                for space in colnames:
                    finale[element].append("HELP")
        
        df_families = pd.DataFrame(finale)
        df_families = df_families.T
        df_families.columns = colnames
       
        return df_families
            
def extractScores(infile, shortNames, normMethod, Z, outputName, familyFile):
    
    #preparing a list of lists with columns for the TFs and each cluster with the scores for them
    scores = {"TF": []}
    transcriptionFactors = {}
    clusterNames = []
    clusterCol = []
    
    if shortNames:
        nameCol = 1
    else:
        nameCol = 0
    
    
    file = open(infile, "r")
    
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
                    

    
    #getting the scores for each TF from each cluster, saving them in the list that was created beforehand

    file.seek(0)   
    for index ,line in enumerate(file):
        #print(currentIndex,line)
        if line.startswith("output_prefix"):
            pass
        else:
            noN = line.split("\n")
            columns = noN[0].split("\t")
            #print(columns)   
                # for element in columns:
                #     try:
                #         columns.remove("")
                #     except:
                #         pass
            
            for i, value in enumerate(clusterCol):
                #print(clusterCol)
                if columns[nameCol] in transcriptionFactors.keys():
                    transcriptionFactors[columns[nameCol]][i] = columns[value]
                else:
                    transcriptionFactors[columns[nameCol]] = []
                    for cluster in clusterNames:
                        transcriptionFactors[columns[nameCol]].append(0)
                    #print(transcriptionFactors)
                    transcriptionFactors[columns[nameCol]][i] = columns[value]
               
    file.close()
                
    for key in transcriptionFactors.keys():
        scores["TF"].append(key)
        
        for index, cluster in enumerate(clusterNames):
            scores[cluster].append(transcriptionFactors[key][index])
    
    
    #print(clusterNames)
    
    #print(scores)    
    #print(transcriptionFactors)
    #checks whether normalization is wanted and which method is supposed to be used
    #normalizes the scores 
    #adjustedScores
    
    possibleNormMethods = ["Sum","None"]
    
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
    
    if familyFile == [""]:
        dfFinal = dfZScores
    else:
        dfFinal = groupFamilies(familyFile, dfZScores)        
    
    colnames = []
    
    for col in dfFinal:
        colnames.append(col)    
    
    dfFinal.to_csv(outputName + '.tsv', sep="\t")    
    
    #CLARION-file generation
    with open(outputName+'.tsv',"r") as csv:
        with open(outputName+'.clarion',"w") as clarion:
            
            clarion.write("!format=Clarion\n#key\tlevel\n#TF\tfeature\n")
            
            for cluster in colnames:
                clarion.write("#"+cluster+"\tsample\n")
        
            for line in csv:
                if line.startswith("\t"):
                    clarion.write("TF"+line)
                else:
                    clarion.write(line)

def main():
    
    inputFile, normMethod, outputName, Z, shortNames, familyFile = cliParser()
    #outfile = open("ClusterComparison.tsv","a")
    print(familyFile)

    extractScores(inputFile[0], shortNames, normMethod, Z, outputName, familyFile[0])


 
if __name__ == '__main__':
    main()