import pandas as pd
import sys

#INPUT: input (input file), method(pearson, kendall, spearman)


filename, input, corr_method, *other = sys.argv

data = pd.read_csv(input, sep='\t')

data.corr(method=corr_method).unstack().sort_values(0, ascending=False).to_csv('{0}_predictions.txt'.format(corr_method), header = None, sep="\t")
