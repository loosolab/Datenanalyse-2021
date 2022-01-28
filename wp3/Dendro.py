import pandas as pd

def main():
    df = pd.read_table("bindetect_distances.txt", sep="\t")
    
    print(df)
    
    families = []
    
    for col in df:
        tmpList = []
        for index in range(len(df.index)):
        
           
            if df.at[index,col] < 0.5:
                tmpList.append(df.iloc[:,index].name)
            
        if tmpList in families:
            pass
        else:
            families.append(tmpList)
    
    #print(tfs)
    
    # families = []
    
    # for key in tfs.keys():
    #      for tf in tfs[key]:
    #          for element in families:
    #              if tf in element:
    #                  pass
    #              else:
    #                  families.append(tfs[key])
    print(families[0],"\n",families[1],"\n",families[2])
        
        
        
    
    
if __name__ == '__main__':
    main()

