#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 00:40:53 2019
@author: ranap

map back interger to string name
"""

import pandas as pd
import networkx as nx
import json
import sys

mapping=json.load(open('mapping.json'))
df= pd.read_csv('coverage_inference.txt',header=-1, delimiter='\t')

map_rev = dict((y,x) for x,y in mapping.items())
df[0] = df[0].map(map_rev)
df.to_csv('mirfluence.output', index=False, header=False)