#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 16:32:24 2021

@author: marlin
"""

import bamnostic as bs

#just for testing how to work with the VMs

#input files, will change to argsparse later on
bam = "/mnt/workspace_stud/stud8/testdata.bam"
tsv = "/mnt/workspace_stud/stud8/clusterIDs.tsv"

#opening the files
allCells = bs.AlignmentFile(bam,"rb")
IDs = open(tsv, "r")

#variables
tabsInIDFile = []
cellIDs = []
clusterIDs = []
clusters = []

#methods (will properly make them later on, just testing how to approach the task)

#splitting the tsv tabs into two lists to have a better ability to work with them
for line in IDs:

        left, right = line.split('\t')

        if right.endswith('\n'):
                right, garbage = right.split('\n')

        if right[0].isdigit():
                right = int(right)

        cellIDs.append(left)
        clusterIDs.append(right)

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
IDs.close()