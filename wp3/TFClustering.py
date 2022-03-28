# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 00:19:06 2022

@author: Moritz Hobein
"""
import argparse
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from os.path import commonprefix

def cliParser():
    parser = argparse.ArgumentParser(description='Clustering transcription factors into families based on motif similarity')
    parser.add_argument('File', nargs=1, help='input file, bindetect_distances.txt')
    args = parser.parse_args()
    
    return args.File

def clusterFamilies(filename):
    data = pd.read_csv(filename, sep="\t")
    
    #clustering, done with scipy
    linkMat = linkage(squareform(data), method="average")
    familyIDs = fcluster(linkMat, 0.5, criterion="distance")	
    
    #grouping the TFs based on their families into a dicitonary
    familyCount = np.unique(familyIDs)
    
    fillFamilies = {}
    
    for ID in familyCount:
        fillFamilies[ID] = []
    
    
    for index, col in enumerate(data):
        fillFamilies[familyIDs[index]].append(col)    
    
    
    familyNames = {}
    
    #naming the clusters/families
    for family in fillFamilies:
        
        for index,TF in enumerate(fillFamilies[family]):
            fillFamilies[family][index] = TF.upper()
        
        #TFs without families
        if len(fillFamilies[family]) == 1:          
            familyNames[fillFamilies[family][0]] = fillFamilies[family]
        
        #TF families where all TFs have the same prefix
        elif len(commonprefix(fillFamilies[family])) > 2 and commonprefix(fillFamilies[family])+" family" not in familyNames.keys(): 
            familyNames[(commonprefix(fillFamilies[family]))+" family"] = fillFamilies[family]
        
        #if there are multiple families with the same prefix (just numbers them)
        elif len(commonprefix(fillFamilies[family])) > 2 and commonprefix(fillFamilies[family])+" family" in familyNames.keys(): 
            counter = 2
            notDone = True
            while(notDone):
                if familyNames.get((commonprefix(fillFamilies[family]))+" family "+str(counter)) == None:
                    familyNames[(commonprefix(fillFamilies[family]))+" family "+str(counter)] = fillFamilies[family] 
                    notDone = False
                else:
                   counter +=1 
        
        #no common name for all TFs in the family possible, take the most common two
        else:    
            tmpNames = []
            etAl = " and others"
            fam = " family"
            
            for index, tester in enumerate(fillFamilies[family]):
                tmpList = []
                for TF in fillFamilies[family]:
                    if len(commonprefix([tester.upper(), TF.upper()])) > 1:
                        tmpList.append(commonprefix([tester.upper(), TF.upper()]))
                       
                tmpNames.append(commonprefix(tmpList))
          
            nameCounter = {}
            for index, TF in enumerate(tmpNames):
                try:
                    shortName, trash = TF.split("_")
                    tmpNames[index] = shortName
                except:
                    pass
                
            for index, TF in enumerate(tmpNames):   
                if TF not in nameCounter:
                    nameCounter[TF] = 1
                else:
                    nameCounter[TF] += 1 
            
            highestCount = 0
            highestName = ""
            secondCount = 0
            secondName = ""
            
            for TF in nameCounter:
               if nameCounter.get(TF) > highestCount:
                   secondName = highestName
                   secondCount = highestCount
                   highestCount = nameCounter.get(TF)
                   highestName = TF
               elif nameCounter.get(TF) > secondCount:
                   secondCount = nameCounter.get(TF)
                   secondName = TF
               else:
                   pass
            
            if len(nameCounter.keys()) > 2:
                name = highestName + ", " + secondName + etAl + fam
            elif secondName == "":
                name = highestName + secondName + fam
            else:
                name = highestName + " and " + secondName + fam
            
            if name not in familyNames.keys():
                familyNames[name] = fillFamilies[family]
            else:
                counter = 2
                notDone = True
                while(notDone):
                    if familyNames.get((name)+" family "+str(counter)) == None:
                        familyNames[name+" family "+str(counter)] = fillFamilies[family] 
                        notDone = False
                    else:
                       counter +=1 
    
    #creating the output file
    #tsv with first column - family name and second column - TFs in the family
    with open("TF_families.tsv","w") as outfile:
        for index, element in enumerate(familyNames):
            if index == 0:
                outfile.write(element + "\t")
            else:
                outfile.write("\n" + element + "\t")
            for i, TF in enumerate(familyNames[element]):
                if i == 0:
                    outfile.write(TF)
                else:
                    outfile.write(","+TF)
                    

def main():
        filename = cliParser()
        clusterFamilies(filename[0])

if __name__ == '__main__':
    main()
