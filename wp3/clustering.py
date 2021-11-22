import bamnostic as bs

#just for testing how to work with the VMs 

#input files, will change to argsparse later on
bam = input(".bam file")
tsv = input(".tsv file for cluster identification")

#opening the files
allCells = bs.AlignmentFile(bam,"rb")
clusterIDs = open(tsv, "r")

for i, read in enumerate(allCells):
  if(i >= 3):
    break
  print(read)

for i, ID in enumerate(clusterIDs):
  if(i >= 3):
    break
  print(ID)

allCells.close()
clusterIDs.close()