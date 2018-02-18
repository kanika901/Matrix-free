% Computes various structural properties of a square matrix.
% Output: nnz, max. no. of non zeros per row etc. 
% Date: February 13, 2018

%A = rand(1000);
%A = load('/Users/kanikas/Documents/research/data_(mat)_matlab_format/mat_files/494_bus.mat');
A = [ 4 1 0; 0 2 1; 0 0 -1];

%1. Dimension 
n = size(A,1);
sprintf('Matrix dimension: %d', n)

%2. No. of non-zeros (nnz)
non_zeros = nnz(A);
sprintf('No. of non-zeros: %d', non_zeros)

nnz_per_row = sum(A ~= 0, 2)
%3. MaxNonzerosPerRow
max_non_zeros_per_row = max(nnz_per_row);
sprintf('Maximum no. of non zeros per row: %d', max_non_zeros_per_row)

%4. MinNonzerosPerRow
min_non_zeros_per_row = min(nnz_per_row);
sprintf('Minimum no. of non zeros per row: %d', min_non_zeros_per_row)

% 5. AvgNonzerosPerRow
avg_non_zeros_per_row = sum(nnz_per_row)/n;
sprintf('Average no. of non zeros per row: %d', avg_non_zeros_per_row)

% DummyRows (No. of rows with only one element)

% 6. AbsoluteNonZeroSum (Sum of the absolute values of all the non zeros in matrix) 
sum_abs_non_zero = sum(sum(abs(A)));
sprintf('Absolute sum of all non zeros: %d', sum_abs_non_zero)

% 7. NumericValueSymmetryV1 (Checks the numerical symmetry: 1 means symmetric, 0 means non-symmetric)
symmetricity = issymmetric(A);
sprintf('Symmetricity: 1 means symmetric, 0 means non-symmetric: %d', symmetricity);

% 8. NonZeroPatternSymmetryV1 (Checks the nonzero pattern symmetry, if symmetric returns 1) 
% if A and A(Transpose) have same nonzero pattern then returns 1)
non_zero_pattern_A = A ~= 0;
non_zero_pattern_A_tran = A' ~= 0;
non_zero_pattern_symmetryV1 = isequal(non_zero_pattern_A, non_zero_pattern_A_tran);
sprintf('Non zero pattern symmetry: %d', non_zero_pattern_symmetryV1)

% 9. Trace
tr = sum(diag(A)); %or trace(A)
sprintf('Trace: %d', tr)

% 10. AbsoluteTrace
abs_trace = sum(diag(abs(A))); 
sprintf('Absolute Trace: %d', abs_trace)

% 11. OneNorm (maximum column sum)
one_norm = max(sum(abs(A)));
sprintf('One Norm: %d', one_norm)

% 12. InfinityNorm (maximum of the absolute row sums)
inf_norm = max(abs(sum(A,2)));
sprintf('Infinity Norm: %d', inf_norm) %or norm(A,'inf')

% 13. FrobeniusNorm
% square root of the sum of the absolute squares of its elements
fro_norm = sqrt(sum(abs(A*A')));
sprintf('Frobenius Norm: %d', fro_norm)

% SymmetricInfinityNorm
% SymmetricFrobeniusNorm
% AntiSymmetricInfinityNorm
% AntiSymmetricFrobeniusNorm
% RowDiagonalDominance 
% ColumnDiagonalDominance 
% RowVariance
% ColumnVariance
% DiagonalAverage
% DiagonalVariance
% DiagonalSign
% DiagonalNonZeros
% lowerBandwidth
% upperBandwidth

