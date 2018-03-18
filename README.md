# Matrix-free approach for matrix feature computation
# Introduction
The traditional way of representing sparse matrices involves storing each nonzero element by using data structures, such as compressed sparse row/column, coordinate, diagonal, or hybrid dense/sparse representations. For certain types of computations, such as nonlinear PDE solution via finite-difference Newton-Krylov methods, where the memory requirements of explicitly storing the sparse matrix exceed available capacity, matrix-free approaches can be used. The Krylov solution of the linearized system is computed by using approximations of matrix-vector products based only on the function computing the current discretized solution approximation at each grid point. Because there is no explicit matrix, it is impossible to compute most of the features used in our ML-based solver selection. Hence, a different set of features must be defined and computed for matrix-free approaches. We present initial results using features based on matrix-free eigenvalue approximation, infinity norm, and structural problem features.
1. The largest eigenvalue is computed by using the Power method (Total iterations = 4).
2. The smallest eigenvalue is computed using Shifted Power method (Total iterations = 4 with double the number of matrix-vector products as in case 1). 
3. These eigenvalues are used for computing the condition number.
 
 # Background
For the Lighthouse project, we have used Anamod and PETSc for computing the matrix features for sparse linear systems. This work encorporates matrix-free computation of features.

# Usage
1. To compute the structural properties of the matrix (Matrix-full):
Script: structural_properties.m
This file writes the properties in a csv file.

2. To compute the eigenvalues with the true eigen values and relative error: 
Script: cond_num_twice_power_mtd_UF.m 
This file writes the properties in a csv file.

The structural and eigenvalues are computed with separate programs. This is because computing eigenvalues are expensive compared to structural values and in the future analysis, if structural properties are enough to make accurate solver decisions, eigenvalues can be skipped altogether.

3. To combine these features into 1 file:
Script: combine_csvs.py
This file combines the structural properties and eigen values based on the matrix name.




