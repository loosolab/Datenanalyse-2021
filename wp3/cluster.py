#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 16:32:24 2021

@authors: moritz, marlin
"""

#import bamnostic as bs
import pysam as ps
import argparse

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
def getSingleClusterBarcodeList (clusterID,combinedIDs):
    cellBarcodesForCluster=[]
    for pair in combinedIDs:
        if pair[1] == clusterID:
            cellBarcodesForCluster.append(pair[0])
    return cellBarcodesForCluster


#Method recives a single cluster ID and a list of CellIDs for this cluster. It the writes a singular BAM file
#containing alls lines that have matching cellIDs
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

#uses tsv file path as input and returns three lists
#cellIDs and clusterIDs contains all IDs in the order they appear, combined IDs is a list of lists with
#two elements containing a pair of IDs that belong together
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
            else:
                continue
            cellIDs.append(left)
            clusterIDs.append(right)
            combinedIDs.append([left,right])
    IDs.close()
    return cellIDs,clusterIDs,combinedIDs
    
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
    cellIDs,clusterIDs,combinedIDs = listifyTSV(args.tsv)
    deduplicatedClusterIDs = deduplicateList(clusterIDs)
    
    #clustering starts here
    for cluster in deduplicatedClusterIDs:
        #getting all cellbarcodes for current cluster
        barcodesForCluster = getSingleClusterBarcodeList(cluster, combinedIDs)
        #writing bam file for current cluster
        writeClusterBam(cluster, barcodesForCluster, args.bam, args.outputDir)
        



if __name__ == "__main__":
    main()