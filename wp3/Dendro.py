import pandas as pd
import numpy as np

def main():
    df = pd.read_table("bindetect_distances.txt", sep="\t")
    
    print(df)
    
    families = {}
    
    filteredCols = []
    filteredScores = []
    
    for i, col in enumerate(df):
        colMeanScore = 0
        tmpList = []
        families[i] =  []
        
        for index in range(len(df.index)):
        
           
            if df.at[index,col] < 0.5:
                tmpList.append(df.iloc[:,index].name)
                colMeanScore = colMeanScore + df.at[index,col]
        
            
        colMeanScore = colMeanScore / len(tmpList)
        
        if tmpList in filteredCols:
            pass
        else:
            filteredCols.append(tmpList)
            filteredScores.append(colMeanScore)
    
    
    #print(filteredScores)
    #print(len(filteredScores))
    #print(df.iloc[:,-2])
    
    
    for index, col in enumerate(df):
        
        #tmpScores = {}
        minMean = 0
        minName = ""
        minPos = 0
        colMeans = []
        indices = []
        scores = []
        
        for i, element in enumerate(filteredCols):
            
            
            if col in element:
                
                scores.append(filteredScores[i])
                indices.append(i)
                        #tmpScores[i] = filteredScores[i]
                        #if filteredScores[i] < minMean:
                            
                           # minimum += filteredScores[i]
                           # minName = col
                           #minPos = index
                           
        tmpMin = 1
        tmpPos = 0
        #print(colMeans)
        for i, score in enumerate(scores):
            if score > 0 and score < tmpMin:
                tmpMin = score
                tmpPos = indices[i]
        
                
            
                       
       # if index <= 10:
            
            #print(tmpScores.keys())
            #print(tmpScores.values())
            #minScore = 1
            #whichmin = ""
            #for index, element in enumerate(tmpScores):
             #   if tmpScores[element] < minScore:
             #      minScore = tmpScores[element]
             #       whichmin = element
            #print(minScore)
            #print(whichmin)    
        
        families[tmpPos].append(col)
                
                
    # for element in families:
    #     for thing in families[element]:
    #         if thing.startswith("TBX"):
                
    #             print(element)
    #print(families[94])
    #print(df.loc[:,"ATF3_MA0605.2"])
        
    #print(tfs)
    
    # families = []
    
    # for key in tfs.keys():
    #      for tf in tfs[key]:
    #          for element in families:
    #              if tf in element:
    #                  pass
    #              else:
    #                  families.append(tfs[key])
        
    outfile = open("TF_familiesv2.txt", "w")
    outfile.write("##FamilyName \t TFs \n")
    
    for index, element in enumerate(families):
        outfile.write(str(index) + "\t")
        
        for TF in families[element]:
            outfile.write(TF+",")
        
        outfile.write("\n")
        
    outfile.close()
    
if __name__ == '__main__':
    main()

