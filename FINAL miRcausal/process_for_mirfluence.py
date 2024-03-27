import pandas as pd
import sys

#INPUT: input (consensus.txt), score (lowest acceptable score), reset_value (what to change scores to for miRfluence)

filename, input, *other = sys.argv
if other:
	try: 
		other[0]=float(other[0])
		score = other[0]
	except ValueError:
		print("Setting score to default value of 0.9")
		score = 0.9
	
	try:
		other[1]=float(other[1])
		reset_value = other[1]
	except ValueError:
		print("Setting reset_value to default value of 0.05")
		reset_value = 0.05
else:
	print("Setting score to default value of 0.9")
	print("Setting reset_value to default value of 0.05")
	score = 0.9
	reset_value = 0.05
list = pd.read_csv(input, sep='\t', header = None) #make sure first row is used for headers
#print(list)
found_division = False
min = 0
max = len(list)

while not found_division:
	
	check_val = int((max+min)/2)
	#print(min, max, score, ":", list[2][check_val], list[2][check_val-1])
	if list[2][check_val]<score and list[2][check_val-1]>=score:
		found_division = True
	elif list[2][check_val] >= score:
		min = check_val
	elif list[2][check_val] < score:
		max = check_val


list = list[:check_val]
list[0] = [gene[1:] for gene in list[0]]
list[1] = [gene[1:] for gene in list[1]]
list[2] = [reset_value] * len(list)
#print (list)
list.to_csv('other_files/filtered_consensus.txt', header = None, index=False, sep="\t")