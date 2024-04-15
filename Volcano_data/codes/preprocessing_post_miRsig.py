import os
import pandas as pd

## Preprocessing before Graph Vizualization ##
'''
path = "D:/scRNA-Seq_data/final_schwartz_networks/"

files = [file for file in os.listdir(path) if file.endswith('network.tsv')]
print(files)

for file in files:
    df = pd.read_csv((path + file), sep="\t", header=None)
    df = df[df[2] >= 0.9]
    df.to_csv((path+file[:-4]+'_final.tsv'), sep='\t', index=False, header=None)
'''

## Preprocessing after Graph Vizualization ##
nodes_path = "D:/scRNA-Seq_data/final_schwartz_networks/"
node_files = [file for file in os.listdir(nodes_path) if (file.startswith("sorted") and file.endswith("90.csv"))]
# print(node_files)

# Keeping only nodes with degree > 5
for file in node_files:
    df = pd.read_csv(nodes_path + file)
    df = df[df['degree'] >= 5]
    df.to_csv(nodes_path + f'{file.split("final")[0]}selected.csv', index=False)

# Read the implant sig files
implants_path = "C:/Users/debnathk/Desktop/study/scRNA-Seq/data/"
implant_files = [file for file in os.listdir(implants_path) if file.endswith("merged.csv")]
# print(implant_files)

# Read the high degree nodes files
hd_node_files = [file for file in os.listdir(nodes_path) if (file.startswith("sorted") and file.endswith("selected.csv"))]
# print(hd_node_files)

# Read the TF file
df_tffile = pd.read_csv('../Rattus_norvegicus_TF.txt', sep='\t')
tf_list = list(df_tffile['Ensembl'])
# print(df_tffile.head())

# Create a union of TF genes and high degree genes
# Genes of Interest
goi_list = []
for file in hd_node_files:
    df_hd_node = pd.read_csv(nodes_path + file)
    goi = set(df_hd_node.iloc[:, 0]).union(set(df_tffile['Ensembl']))
    goi_list.append(list(goi))

# Save the implant deg files with goi (hd genes + TF genes)
for (file, goi) in zip(implant_files, goi_list):
    df_implants = pd.read_csv(implants_path + file)
    df_implants = df_implants[df_implants['ENSEMBL ID'].isin(goi)]
    df_implants['ENSEMBL ID'].to_csv(implants_path + f'{file.split(".csv")[0]}_goi_genes.csv', index=False)
    df_implants.to_csv(implants_path + f'{file.split(".csv")[0]}_goi.csv', index=False)

