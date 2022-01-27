import pandas as pd

df = pd.read_table("bindetect_distances.txt", sep="\t")

print(df)

tfs = {}

for col in df:
	tfs[col] = []

	for index in range(len(df.index)):
		if df.at[index,col] < 0.4:
			tfs[col].append(col)

print(tfs)
