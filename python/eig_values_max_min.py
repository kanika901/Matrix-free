#!/usr/bin/env python3 

'''
Date: February 8, 2018 
Author: Kanika Sood

Computes the largest and smallest eigen value using eigvals() [absolute values]
Only supports square matrices for now.
To run: python3 eig_values_max_min.py 
'''

import numpy as np
from scipy.sparse.linalg import eigs
from scipy.io import loadmat
from scipy.sparse import identity
from os.path import basename

class EigenValues:

	def changeMatrixFormat(self, B):
		A = loadmat(B)
		matrix = A['Problem']['A']
		A = np.array(matrix, dtype = object)
		A = A[0][0].toarray() #converting scipy csc.csc matrix to numpy ndarray
		self.computeEigenValues(A)
		return A

	def trueEigenValues(self, A):
		print('Computing true eigen values.....')
		lamda_max = max(abs(np.linalg.eigvals(A)))
		lamda_min = min(np.linalg.eigvals(A))
		#print('Eigen values: ', np.linalg.eigvals(A))
		return lamda_max, lamda_min

	def largestEigenValue(self, A):
		#using power method for largest eigen value with MV products and fixed no. of iterations (4)
		print('Using power method to compute largest eigen value....')
		x, y = A.shape
		if x != y:
			print('Sorry, this program only supports square matrices for now')
		elif x == y:
			n = x
		w = np.random.rand(n,1)
		w = w/np.linalg.norm(w)
		cntr = 0
		itr = 0
		while itr < 4: 
			cntr += 1
			itr += 1
			b = A * w #matrix-vector product --> cpk
			lamda = np.max(abs(b)) #max(abs(b))
			lamda_with_sign = np.amax(b)
			if lamda == lamda_with_sign:
				lamda_sign = 1
			else: 
				lamda_sign = -1
			lamda = lamda * lamda_sign
			w = b/lamda
		return lamda, n

	def smallestEigenValue(self, A, lamda_max_A, n): 
		#using shifted power method for smallest eigen value with MV products and fixed no. of iterations (4)
		print('Using shifted power method to compute smallest eigen value....')
		I = np.identity(n)
		B = A - (lamda_max_A * I)
		lamda_max_B, n = self.largestEigenValue(B)
		lamda_min_A = (lamda_max_B) + (lamda_max_A)
		return lamda_min_A

	def computeEigenValues(self, A):
		lamda_max_A, n = self.largestEigenValue(A)
		print('Largest eigen value from power method: ', lamda_max_A)
		lamda_min_A = self.smallestEigenValue(A, lamda_max_A, n)
		print('Smallest eigen value using shifted power method: ', lamda_min_A)
		#self.conditionNum(lamda_max_A, lamda_min_A, A)
		lamda_max_true, lamda_min_true = self.trueEigenValues(A)
		print('True largest eigen value: ', lamda_max_true)
		print('True smallest eigen value: ', lamda_min_true)


if __name__ == '__main__':

	eigvals = EigenValues()
	#A = np.array([[21, 11, 7], [ 1, 18, 10], [4, 8, 11]])
	#A = '/Users/kanikas/Documents/research/data_(mat)_matlab_format/mat_files/494_bus.mat'
	#A = '/Users/kanikas/Documents/research/data_(mat)_matlab_format/mat_files/2cubes_sphere.mat'
	#A = '/Users/kanikas/Documents/research/data_(mat)_matlab_format/mat_files/m3plates.mat'
	A = '/Users/kanikas/Documents/research/data_(mat)_matlab_format/mat_files/barth5.mat'
	print('Computing eigen values for matrix: ', basename(A))
	matrix = eigvals.changeMatrixFormat(A)

	

