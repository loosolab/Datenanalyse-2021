# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 12:02:36 2022

@author: MoritzHobein
"""

import pandas as pd
import matplotlib.pyplot as plt
import argparse

def cliParser():
    parser = argparse.ArgumentParser(description='Generates a plot and a file, where and how often the most important transcription factors per cluster appear')
    parser.add_argument('File', nargs="*", help='output of DefiningTF.py (tsv)')
    parser.add_argument('-f', '--families', dest='Families', default="idontwantanyfamilies", help='output of TFClustering.py (tsv)')
    args = parser.parse_args()
    
    return args.File, args.Families

def extractData(inputFiles, familyFile):
    #inputFiles = ["testdaten/liverDTF.txt","testdaten/lungDTF.txt"]
    
    colorMapping = {}
    TFs = {}
    
    for file in inputFiles:
        df = pd.read_csv(file, sep="\t")
    
        for col in df:
            factors = df[col].tolist()
            colorMapping[col]=factors
            
            for tf in factors:
                if tf in TFs.keys():
                    pass
                else:
                    TFs[tf] = []
        
        for col in df:
            factors = df[col].tolist()
            
            for tf in factors:
                TFs[tf].append(col)
                
                
    return colorMapping, TFs


def bubbleSize(TFs, familyFile):
    familySizes = {}
    
    with open(familyFile, "r") as file:
        for line in file:
            familyName, members = line.split("\t")
            members = members.split(",")
            memberCount = len(members)
            
            if familyName in TFs.keys():
                familySizes[familyName] = memberCount*3 + 15
    
    return familySizes 

#TODO überschrift, label, etc
def generatePlot(colorMapping, TFs, familySizes): 
    
    tfx = []
    clustery = []
    colorsz = []
    bubbleSizes = []
    
    for index, tf in enumerate(TFs):
        
        for cluster in TFs[tf]:
            if type(tf) == int:
                pass
            else:
                colorValue = len(TFs[tf])
                tfx.append(tf)
                clustery.append(cluster)
                colorsz.append(colorValue)
                bubbleSizes.append(familySizes[tf])
                
    
    plt.scatter(tfx, clustery, c=colorsz, cmap = "Set1", s=bubbleSizes)
    
    plt.xticks(rotation=30,ha='right', size=5)
    plt.colorbar()
    plt.show()
 
    
def createTSV(TFs):
       
    with open("similarity.tsv","w") as outfile:
        
        outfile.write("TF\tcelltypes\n")
        
        for tf in TFs:
            outfile.write(str(tf) + "\t")
            
            for cluster in TFs[tf]:
                outfile.write(str(cluster) + "\t")
            outfile.write("\n")
    
                
        
def main():
    
    inputFiles, familyFile = cliParser()
    
    colorMapping, TFs = extractData(inputFiles, familyFile)
    
    if familyFile == "idontwantanyfamilies":
        
        familySizes = {}
        
        for tf in TFs:
            familySizes[tf] = 10
    
    else:
        familySizes = bubbleSize(TFs, familyFile)
    
    generatePlot(colorMapping, TFs, familySizes)
    
    createTSV(TFs)


if __name__ == '__main__':
    main()
      
    

