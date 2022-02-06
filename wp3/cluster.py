#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 16:32:24 2021

@authors: moritz, marlin
"""

import pysam as ps
import argparse

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
#no longer in use may be removed
def getSingleClusterBarcodeList (clusterID,combinedIDs):
    cellBarcodesForCluster=[]
    for pair in combinedIDs:
        if pair[1] == clusterID:
            cellBarcodesForCluster.append(pair[0])
    return cellBarcodesForCluster


#Method recives a single cluster ID and a list of CellIDs for this cluster. It the writes a singular BAM file
#containing alls lines that have matching cellIDs
#No longer in use may be removed
def writeClusterBam (clusterID,cellIDsForCluster,sourceFilePath,outputDir):
    #debug print
    print("working on cluster: "+str(clusterID))
    #print("Debug, cellbacodes for this cluster: "+str(cellIDsForCluster))
    sourceFile = ps.AlignmentFile(sourceFilePath,"rb")
    clusterFile = ps.AlignmentFile(outputDir+"cluster"+str(clusterID)+".bam","wb",template=sourceFile)
    for read in sourceFile.fetch():
        if read.has_tag('CB'):
            if read.get_tag('CB') in cellIDsForCluster:
                #print("Debug, matching read for tag: "+str(read.get_tag('CB')))
                clusterFile.write(read)
    sourceFile.close()
    clusterFile.close()

#hopefully faster way to write single cluster bam files
def writeClusters (deduplicatedClusterIDs,cellClusterDict,sourceFilePath,outputDir):
    #creating files
    sourceFile = ps.AlignmentFile(sourceFilePath,"rb")
    clusterFileDict = {}
    #creating new files
    for cluster in deduplicatedClusterIDs:
        singleClusterFile = ps.AlignmentFile(outputDir+"cluster"+str(cluster)+".bam","wb",template=sourceFile)
        clusterFileDict.update({cluster:singleClusterFile})
    #scanning source for lines and writing to single cluster files
    for read in sourceFile.fetch():
        if read.has_tag('CB:Z'):
            try:
                clusterID = cellClusterDict[read.get_tag('CB:Z')]
            except:
                #uncomment next line if you wish to be made aware of Reads that don't have a cluster assigned to them.
                #print("Line with CB:Z: "+read.get_tag('CB:Z')+" does not have a cluster assignment and will be ignored.")
                continue
            clusterFileDict[clusterID].write(read)
    for file in clusterFileDict.values():
        file.close()
    sourceFile.close()

#uses tsv file path as input and returns three lists
#cellIDs and clusterIDs contains all IDs in the order they appear, combined IDs is a list of lists with
#two elements containing a pair of IDs that belong together
#No longer in use may be removed
def listifyTSV(tsvPath):
    cellIDs = []
    clusterIDs = []
    combinedIDs = []
    cellClusterDict = {} #{CellID:ClusterID}
    IDs = open(tsvPath, "r")
    for line in IDs:

            left, right = line.split('\t')

            if right.endswith('\n'):
                    right, garbage = right.split('\n')

            if right[0].isdigit():
                    right = int(right)
            else:
                continue
            cellIDs.append(left)
            clusterIDs.append(right)
            combinedIDs.append([left,right])
            cellClusterDict.update({left:right})
    IDs.close()
    return cellIDs,clusterIDs,combinedIDs,cellClusterDict
    
#returns a dictionary that contains {Cellbarcode:ClusterID}
def tsvToDict(tsvPath):
    assignmentDict = {}
    clusterList = []
    tsv = open(tsvPath, "r")
    for line in tsv:
        barcode, clusterID = line.split("\t")
        assignmentDict.update({barcode.strip():clusterID.strip()})
        clusterList.append(clusterID.strip())
    tsv.close()
    return assignmentDict, clusterList


#generates a text file that can be pasted into the TOBIAS snakemake pipline config file
def generateSnakemakeInput(outputDir, clusterIDs):
    f = open(outputDir+"snakemakeIn.txt","a")
    for ID in clusterIDs:
        f.write(" cluster"+str(ID)+": ["+outputDir+"cluster"+str(ID)+".bam"+"]\n")
    f.close()
    
def main():
    #Argparser
    aParser = argparse.ArgumentParser(description='This is a tool for splitting a .bam file into multiple .bam files.')
    aParser.add_argument('-Bam', '-b',dest='bam', nargs='?', help='Use this to set the input Path for the .bam file.',
                         default="inputWP3/testdata.bam")
    aParser.add_argument('-Tsv', '-t',dest='tsv', nargs='?', help='Use this to set the input Path for the .tsv file,'+
                         'containing the cluster Assingnments and cell barcodes.',
                         default="inputWP3/clusterIDs.tsv")
    aParser.add_argument('-Out', '-o',dest='outputDir', nargs='?', help='Use this to set the Path to a directory,'+
                         'that should be used to save this programs output.',
                         default="outputWP3/")
    args=aParser.parse_args()
    #read tsv file and convert it into lists
    #cellIDs,clusterIDs,combinedIDs,cellClusterDict = listifyTSV(args.tsv)

    clusterDict, clusterIDs = tsvToDict(args.tsv)
    deduplicatedClusterIDs = deduplicateList(clusterIDs)
    
    writeClusters(deduplicatedClusterIDs, clusterDict, args.bam, args.outputDir)
    
    #clustering starts here
    # for cluster in deduplicatedClusterIDs:
    #     #getting all cellbarcodes for current cluster
    #     barcodesForCluster = getSingleClusterBarcodeList(cluster, combinedIDs)
    #     #writing bam file for current cluster
    #     writeClusterBam(cluster, barcodesForCluster, args.bam, args.outputDir)
    
        
    #generating file for TOBIAS Snakemake pipeline
    generateSnakemakeInput(args.outputDir, deduplicatedClusterIDs)

if __name__ == "__main__":
    main()