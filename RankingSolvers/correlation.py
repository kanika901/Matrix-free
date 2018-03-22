#!/usr/bin/env python
'''
Created on March 21, 2018
Usage: python3 correlation.py
Note: The csv input file should have class labels as integer and not strings
'''

import pandas as pd
import numpy as np 
from scipy.stats.stats import pearsonr, spearmanr
import itertools


file_location = '/Users/kanikas/Documents/MatrixFree/Matrix-free/MLdata/filesUsedForResults/'
filename = 'correlation_RS1_mfree_UF_arya_p72_final_45.csv'
data = pd.read_csv(file_location+filename)

correlations_pearson = {}
correlations_spearman = {}

columns = data.columns.tolist()
print(columns)

for col1, col2 in itertools.combinations(columns,2):
	correlations_pearson[col1 + ' <-> ' + col2] = pearsonr(data.loc[:,col1], data.loc[:,col2])
	correlations_spearman[col1 + ' <-> ' + col2] = spearmanr(data.loc[:,col1], data.loc[:,col2])

result_pearson = pd.DataFrame.from_dict(correlations_pearson, orient='index')
result_spearman = pd.DataFrame.from_dict(correlations_spearman, orient='index')
result_pearson.columns = ['PCC', 'p-value']
result_spearman.columns = ['PCC', 'p-value']

print('-----Pearson correlation-----')
print(result_pearson.sort_index())

result_spearman.sort_index().to_csv('correlations_spearman.csv')

print('-----Spearman correlation-----')
print(result_spearman.sort_index())
result_pearson.sort_index().to_csv('correlations_pearson.csv')




