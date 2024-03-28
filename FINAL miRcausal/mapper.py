#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 23:57:05 2019

@author: ranap
This code convert a nodename of a network to integer


"""

import pandas as pd
import networkx as nx
import json
import sys

df = pd.read_csv(str(sys.argv[1]), header=-1 )
G = nx.from_pandas_edgelist(df,0,1,2)

mapping = dict(zip(G.nodes(), range(1, G.number_of_nodes()+1)))
G1 = nx.relabel_nodes(G, mapping)
nx.write_edgelist(G1, "test.edgelist", data=False, delimiter='\t')

df=pd.read_csv('test.edgelist', delimiter='\t', header=-1)
print(df.shape)
df[2] = 0.25
df.to_csv('test.edgelist', sep='\t', index=False, header=False)  

with open('mapping.json', 'w') as fp:
    json.dump(mapping, fp)

