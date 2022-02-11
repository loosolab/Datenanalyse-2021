# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 00:19:06 2022

@author: Moritz Hobein
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.cluster.hierarchy as sch
from sklearn.datasets import load_iris
from sklearn.cluster import AgglomerativeClustering

data = pd.read_csv("bindetect_distances.txt", sep="\t")

df = data.iloc[0:10].copy()

print(df)

dendrogram = sch.dendrogram(sch.linkage(df, method  = "ward"))
plt.title('Dendrogram')
plt.xlabel('TFs')
plt.ylabel('Euclidean distances')
plt.show()

# agg = AgglomerativeClustering(n_clusters=3,affinity = 'euclidean', linkage = 'ward',
#                               #distance_threshold=0.5
#                               )
                              
# labels = agg.fit_predict(df)

print(sch.linkage(df,method = "ward"))


