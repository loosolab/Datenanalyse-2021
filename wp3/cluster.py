#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 16:32:24 2021

@author: marlin
"""

import bamnostic as bs



#methods (will properly make them later on, just testing how to approach the task)

#for deduplicating lists, will be used to calculate number of distinct clusters and output files
def deduplicateList (listWithDuplicates):
    deduplicated=[]
    for entry in listWithDuplicates:
        if entry in deduplicated:
            continue
        else:
            deduplicated.append(entry)
    return deduplicated

#for getting a list of all cell Barcodes, that belong to the provided cluster.
#TODO still needs an input of cellbarcodes
def getSingleClusterBarcodeList (clusterID):
    cellBarcodesForCluster=[]
    return cellBarcodesForCluster


#Method recives a single cluster ID and a list of CellIDs for this cluster. It the writes a singular BAM file
#containing alls lines that have matching cellIDs
def writeClusterBam (clusterID,cellIDsForCluster,sourceFilePath,outputDir):
    sourceFile = bs.AlignmentFile(sourceFilePath,"rb")
    clusterFile = bs.AlignmentFile(outputDir+"cluster"+str(clusterID)+".bam","wb",template=sourceFile)
    for read in sourceFile.fetch():
        #TODO not sure if refID is the right value to call needs testing
        if read.refID in cellIDsForCluster:
            clusterFile.write()
    sourceFile.close()
    clusterFile.close()

#uses tsv file path as input and returns three lists
def listifyTSV(tsvPath):
    cellIDs = []
    clusterIDs = []
    combinedIDs=[]
    IDs = open(tsvPath, "r")
    for line in IDs:

            left, right = line.split('\t')

            if right.endswith('\n'):
                    right, garbage = right.split('\n')

            if right[0].isdigit():
                    right = int(right)

            cellIDs.append(left)
            clusterIDs.append(right)
            combinedIDs.append([left,right])
    IDs.close()
    return cellIDs,clusterIDs,combinedIDs
    
def main():
    #just for testing how to work with the VMs

    #input files, will change to argsparse later on
    bam = "inputWP3/testdata.bam"
    tsv = "inputWP3/clusterIDs.tsv"
    outputDir="outputWP3/"

    #opening the files
    allCells = bs.AlignmentFile(bam,"rb")
    
    #reading tsv
    cellIDs,clusterIDs,combinedIDs = listifyTSV(tsv)

    #variables
    tabsInIDFile = []
    clusters = []
    
    #splitting the tsv tabs into two lists to have a better ability to work with them


    #putting both tabs into one list
    tabsInIDFile.append(cellIDs)
    tabsInIDFile.append(clusterIDs)
    print(tabsInIDFile)

    cellIDs.pop(0)
    clusterIDs.pop(0)
    print(max(clusterIDs))

    for i in range(1,max(clusterIDs)+1):
            clusters.append(i)

    print(clusters)

    #for i, read in enumerate(allCells):
     # if(i >= 3):
      #  break
      #print(read)

    #for i, ID in enumerate(clusterIDs):
     # if(i >= 3):
      #  break
     # print(ID)

    allCells.close()


if __name__ == "__main__":
    main()