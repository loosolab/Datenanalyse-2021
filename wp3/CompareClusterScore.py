# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 20:09:09 2021

@author: Moritz Hobein
"""

import argparse
import pandas as pd

#parser for using the tool via command-line interface
def cliParser():
    parser = argparse.ArgumentParser(description='Graphic Cluster Footprint Score Comparison')
    parser.add_argument('Files', nargs="*", help='input file, bindetect_results.txt. If supplying the script with multiple files, rename the cluster names according to tissue to avoid conflicts (i.e. with sed)')
    parser.add_argument('-n', '--normalize', dest='Norm', default='Sum', help='Method of normalization of the data.')
    parser.add_argument('-o', '--outputName', dest='Custom_filename', default="ClusterComparison", help='Sets a custom name for saving the output files')
    parser.add_argument('-z', '--ZScore', dest='Z', default=True, action='store_false', help='if used, skips Z Scores calculation for each TF')
    parser.add_argument('-s', '--shortNames', dest='shortNames', default=False, action='store_true', help='If used, saves the short Name of each TF instead of the output_prefix. Does not work with family clusters')
    parser.add_argument('-f', '--families', dest='families', nargs=1, default="idontwantanyfamilies", help='Family cluster file to group TF families together')
    args = parser.parse_args()
    
    return args.Files, args.Norm, args.Custom_filename, args.Z, args.shortNames, args.families

#normalizing the scores to the same sum
def normalizeToSum(scores, fileNames):
     
    normalizationFactors = []
    index = 0
    
    for element in scores:
        if element != "TF":
            normalizationFactors.append(0)
            
            for value in scores[element]:
                normalizationFactors[index] = normalizationFactors[index] + float(value)
            index = index + 1  
    
    df = pd.DataFrame(scores)
    df.set_index('TF', inplace=True)
    
    df = df.astype(float)
    
    for i, element in enumerate(fileNames):
        df[element] =  df[element]*(100/normalizationFactors[i])
    
    
    return df
       
#calculating z scores from the normalized footprinting scores
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

    
    dfZ = dfZ_transposed.transpose()

    return dfZ

#compressing the TFs into clusters    
def groupFamilies(familyFile, df):
   
    with open(familyFile[0], "r") as file:
        
        colnames = []
        
        for col in df:
            colnames.append(col)
        
        dfT = df.transpose()
        
        families = {}

        for line in file:
            familyName, familyMembers = line.split("\t")
            familyMembers= familyMembers.split("\n")
            familyMembers = familyMembers[0].split(",")
            
            for index, member in enumerate(familyMembers):
                familyMembers[index] = member.upper()
            
            families[familyName] = familyMembers
        

        familyScores = {}
        
        for familyName in families.keys():
            tmpScores = []
            memberCounter = 0
            
            
            for col in dfT:
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
         
            meanScores = []
            
            for score in tmpScores:
                meanScores.append(score/memberCounter)
            
            familyScores[familyName] = meanScores
        
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

#extracting the scores from the input files and creating output files containing them           
def extractScores(infile, shortNames, normMethod, Z, outputName, familyFile):
    
    #preparing a list of lists with columns for the TFs and each cluster with the scores for them
    scores = {"TF": []}
    transcriptionFactors = {}
    clusterNames = []
    
    #using full names or short names (family clustering works with short names only)
    if shortNames:
        nameCol = 1
    else:
        nameCol = 0
    
    #for each file, get the relevant columns for data extraction
    for tissue in infile:
        clusterCol = []
    
        file = open(tissue, "r")
        
        for index, line in enumerate(file):

            if index == 0:
                noN, garbage = line.split("\n")
                elements = noN.split("\t")

                for i, col in enumerate(elements):
                    if col.endswith("_mean_score"):
                        name = col.split("_mean_score")
                        clusterNames.append(name[0])
                        clusterCol.append(i)
                        scores[name[0]] = []
                        
    

        #getting the scores for each TF from each cluster, saving them in the list that was created beforehand
    
        file.seek(0)   
        for index ,line in enumerate(file):
            if line.startswith("output_prefix"):
                pass
            else:
                noN = line.split("\n")
                columns = noN[0].split("\t")
                
                #writing down the scores for each TF and cluster
                for i, value in enumerate(clusterCol):

                    if columns[nameCol] in transcriptionFactors.keys():
                        
                        #assigning the score of the first cluster to a TF
                        if transcriptionFactors[columns[nameCol]][i] == "EMPTY":
                            transcriptionFactors[columns[nameCol]][i] = columns[value]
                        #adding the scores of other clusters
                        else:
                            transcriptionFactors[columns[nameCol]].append(columns[value])
                    
                    #if the TF has not been visited yet, mke it possible to assign a score
                    else:
                        transcriptionFactors[columns[nameCol]] = []
                        for cluster in clusterNames:
                            transcriptionFactors[columns[nameCol]].append("EMPTY")

                        transcriptionFactors[columns[nameCol]][i] = columns[value]          
        file.close()
    
    #creating a dictionary containing the scores for each TF and cluster to create a DataFrame later
    for key in transcriptionFactors.keys():
        scores["TF"].append(key)
        
        for index, cluster in enumerate(clusterNames):
            scores[cluster].append(transcriptionFactors[key][index])
    
    
    #normalizing the data 
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

    #calculating the z score
    if Z:        
        dfZScores = toZScore(dfNormScores, clusterNames)  
    else:
        dfZScores = dfNormScores
    
    #grouping TFs into families
    if familyFile == "idontwantanyfamilies":
        dfFinal = dfZScores
    else:
        dfFinal = groupFamilies(familyFile, dfZScores)        
    
    #necessary for CLARION file generation
    colnames = []
    
    for col in dfFinal:
        colnames.append(col)    
    
    #tsv file generation (necessary to create the CLARION file, also useful if the CLARION header is not wanted)
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

    extractScores(inputFile, shortNames, normMethod, Z, outputName, familyFile)


 
if __name__ == '__main__':
    main()
