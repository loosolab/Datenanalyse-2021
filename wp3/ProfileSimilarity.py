# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 12:02:36 2022

@author: MoritzHobein
"""

import pandas as pd
import matplotlib.pyplot as plt

inputFiles = ["testdaten/liverImp.tsv"]

TFs = {}

for file in inputFiles:
    df = pd.read_csv(file, sep="\t")

    for col in df:
        factors = df[col].tolist()
        
        for tf in factors:
            if tf in TFs.keys():
                pass
            else:
                TFs[tf] = []
    
    for col in df:
        factors = df[col].tolist()
        
        for tf in factors:
            TFs[tf].append(col)


for tf in TFs:
    for cluster in TFs[tf]:

        plt.scatter(tf, cluster)

plt.legend()
plt.show()
    
with open("similarity.tsv","w") as outfile:
    
    outfile.write("#TF\tcelltypes\n")
    
    for tf in TFs:
        outfile.write(str(tf) + "\t")
        
        for cluster in TFs[tf]:
            outfile.write(str(cluster) + "\t")
        outfile.write("\n")
    
                
        
        
    

