import pandas as pd
import sys
import pickle
#INPUT: input (miRfluence output), percent (percent of top to keep)

filename, input, *other = sys.argv

if other:
	try:
		other[0] = float(other[0])
		percent = other[0]
	except ValueError:
		print("Setting percent to default value of 0.1")
		percent = 0.1
else:
	print("Setting percent to default value of 0.1")
	percent = 0.1

list = pd.read_csv(input, sep='\t', header = None) #make sure first row is used for headers
pickle.dump(list, open("other_files/list.pkl", 'wb'))

list = list.sort_values(1, ascending=False)
#print("sortedlist", list)
col = list[0][:round(len(list)*percent)]
#print("00000\n\n" , list[0])
#print("00000\n\n" , list[0][:round(len(list)*percent)])
#print (round(len(list)*percent))

conversions = pickle.load(open("other_files/gene_conversions.pkl", "rb"))
#print (conversions)
new_col = []
for val in col:
	new_col.append(conversions[val])
#print(new_col)
pd.DataFrame(new_col).to_csv('influential_mirna.txt', header = None, index=False, sep="\t")
