import pandas as pd
import pickle
import sys

# map genename to string. Save the genename to string mapping


filename, input, *other = sys.argv

data = pd.read_csv(input, sep='\t') #make sure first row is used for headers

headers_to_remove = []

data = data[[header for header in data if data[header].nunique() > 2 ]]

conversions = {key:val for key, val in enumerate(data)}

pickle.dump(conversions, open("other_files/gene_conversions.pkl", 'wb'))
data.to_csv('other_files/cleaned_matrix.txt', header = ["G"+str(i) for i in conversions.keys()], index=False, sep="\t")