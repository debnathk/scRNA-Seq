#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 00:40:53 2019
@author: ranap

map back interger to string name
"""

import pandas as pd
import pickle

df = pd.read_csv("mirsig_network.txt", sep='\t', header=None)
file = open("other_files/gene_conversions.pkl",'rb')
mapping = pickle.load(file)
file.close()
mapping = {'G' + str(k): v for k, v in mapping.items()}

df[0] = df[0].map(mapping)
df[1] = df[1].map(mapping)

df.to_csv("mirsig_network_final.tsv", index=False, header=False, sep='\t')
