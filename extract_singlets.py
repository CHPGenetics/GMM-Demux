import pandas as pd

from sys import argv

df = pd.read_csv(argv[1], index_col=0)

index_list = df.loc[(df['Cluster_id'] > 0) & (df['Cluster_id'] < 5)].index.tolist()

with open('SSD.csv', 'w') as f:
    for item in index_list:
        f.write("%s\n" % item)

print(df.loc[index_list])
