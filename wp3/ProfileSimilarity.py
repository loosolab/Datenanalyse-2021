# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 12:02:36 2022

@author: MoritzHobein
"""

import pandas as pd
import matplotlib.pyplot as plt

<<<<<<< HEAD
inputFiles = ["testdaten/liverImp.tsv"]
=======

inputFiles = ["testdaten/liverDTF.txt","testdaten/lungDTF.txt"]
>>>>>>> d8881417cbb65c42ff5dd11be4f855a51ddc296a

TFs = {}

for file in inputFiles:
    df = pd.read_csv(file, sep="\t")
<<<<<<< HEAD

=======
    #print(df)
    
>>>>>>> d8881417cbb65c42ff5dd11be4f855a51ddc296a
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
<<<<<<< HEAD

=======
            
>>>>>>> d8881417cbb65c42ff5dd11be4f855a51ddc296a

for tf in TFs:
    for cluster in TFs[tf]:

        plt.scatter(tf, cluster)

plt.legend()
plt.show()
<<<<<<< HEAD
=======

>>>>>>> d8881417cbb65c42ff5dd11be4f855a51ddc296a
    
with open("similarity.tsv","w") as outfile:
    
    outfile.write("TF\tcelltypes\n")
    
    for tf in TFs:
        outfile.write(str(tf) + "\t")
        
        for cluster in TFs[tf]:
            outfile.write(str(cluster) + "\t")
        outfile.write("\n")
    
                
        
        
    

