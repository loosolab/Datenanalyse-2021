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
    parser.add_argument('-f', '--families', dest='Families', default="idontwantanyfamilies", help='output of TFClustering.py (tsv). Required to match bubble size to family size')
    parser.add_argument('-o', '--output', dest='Output', default="similarity.tsv", help='output of TFClustering.py (tsv)')
    parser.add_argument('-t', '--title', dest='Title', default="Bubble Plot of the Defining Transcription Factors per Cluster", help='setting a custom plot title')
    args = parser.parse_args()
    
    return args.File, args.Families, args.Output, args.Title

#takes the data provided by the input files and puts them into dictionaries to generate the plot with
def extractData(inputFiles, familyFile):
    
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

#sets the bubble size of each family in the plot proportional to the family size
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

#plots the data
#tfs on x axis, clusters on y axis, colors according to overlaps between cluster profiles, bubble size according to family size
def generatePlot(colorMapping, TFs, familySizes, title, outputName): 
    
    if outputName == "similarity.tsv":
        plotFileName = "SimilarityBubblePlot.png"
    else:
        withoutFiletype, Filetype = outputName.split(".")
        plotFileName = "withoutFiletype"+"BubblePlot.png"
    
    tfx = []
    clustery = []
    colorsz = []
    bubbleSizes = []
    
    #filling out lists with data to use them to generate the plot afterwards
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
    
    plt.xlabel("Transcription Factors")
    plt.ylabel("Clusters")
    plt.title(title)
    plt.xticks(rotation=30,ha='right', size=5)
    plt.colorbar()
    
    fig = plt.gcf()
    
    fig.savefig(plotFileName)
    plt.show()
 
#creates a file that contains each TF in the left column and all clusters that show this TF in their top TFs in the right column    
def createTSV(TFs, outputName):
       
    with open(outputName,"w") as outfile:
        
        outfile.write("TF\tcelltypes\n")
        
        for tf in TFs:
            outfile.write(str(tf) + "\t")
            
            for cluster in TFs[tf]:
                outfile.write(str(cluster) + "\t")
            outfile.write("\n")
    
                
        
def main():
    
    inputFiles, familyFile, outputName, title = cliParser()
    
    colorMapping, TFs = extractData(inputFiles, familyFile)
    
    #if there is no family file, creates a dictionary that sets all bubbles to the same size
    if familyFile == "idontwantanyfamilies":
        
        familySizes = {}
        
        for tf in TFs:
            familySizes[tf] = 10
    
    else:
        familySizes = bubbleSize(TFs, familyFile)
    
    generatePlot(colorMapping, TFs, familySizes, title, outputName)
    
    createTSV(TFs, outputName)


if __name__ == '__main__':
    main()
      
    

