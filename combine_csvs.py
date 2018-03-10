#!/usr/bin/env python3

'''
Date: March 9, 2018
Author: Kanika Sood
Takes 2 csv files and merges them on the basis of the matrix name. 
Note: The matrix name should be the first column of both the input csv files
'''

import csv 

filename1 = '/Users/kanikas/Documents/MatrixFree/Matrix-free/StucturalProperties/struct_properties_all_matr_final_v3.csv'
filename2 = '/Users/kanikas/Documents/MatrixFree/Matrix-free/ConditionNo/eig_vals_cond_num_iter4_rel_err_all_matr_final_v2.csv'

outfile = open('/Users/kanikas/Documents/MatrixFree/Matrix-free/combined_structural_eigen.csv', 'w+')
out_writer = csv.writer(outfile)

infile1 = open(filename1, 'r')
infile2 = open(filename2, 'r')

header1 = next(infile1).split(',')
header2 = next(infile2).split(',')

out_writer.writerow(header1 + header2)

for row1 in infile1:
		row1 = row1.split(',')
		for row2 in infile2:
				row2 = row2.split(',')
				if str(row1[0]) == str(row2[0]): #combine csv based on matrix name
					out_writer.writerow(row1 + row2)
					print('Matrix: ', row1[0])

		infile2 = open(filename2, 'r')
		header2 = next(infile2)

infile2.close()
infile1.close()
outfile.close()