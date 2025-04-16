import pandas as pd
import sys
print("---------Merging files-----------")
files=list(sys.argv[1:])
for i,file in enumerate(files):
    f=pd.read_csv(file)
    if i==0:
        df=f
    else:
        df=pd.concat([df,f])
df.to_csv("merged_A01_A10_matchedonly.csv")

print("---------Finished-----------")