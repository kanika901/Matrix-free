% Computes RS1 structural properties of a square matrix.
% test set for ML classification
% Output: nnz, max. no. of non zeros per row etc. 
% Date: March 21, 2018

cd '/Volumes/Kank/MatrixFree'
file_location = '/Volumes/Kank/MatrixFree/TestMatrices';
%out_file = fopen('/Volumes/Kank/MatrixFree/RS1test_set2.csv','a') ;
%cd '/Users/kanikas/Documents/MatrixFree/Matrix-free'
%file_location = '/Users/kanikas/Documents/MatrixFree/Matrix-free/TestMatrices';

all_files = dir(file_location);
all_names = { all_files.name };
global in_file;
out_file = fopen('/Volumes/Kank/MatrixFree/RS1test_set2.csv','a') ;
fprintf(out_file,'Matrix name, MinNonzerosPerRow, NonZeroPatternSymmetryV1, InfinityNorm, ColumnVariance, ColumnVarianceNonZeros, lowerBandwidth, DiagonalNonZeros, DiagonalAverage \n');

for j = 4: length(all_names)
    in_file = strcat(file_location,'/', all_files(j).name);
    [path, mat_name, ext] = fileparts(in_file);
    A = load(in_file);
    matname = strsplit(mat_name, ',');
    mat_name = matname{:};
    A = getfield(A,mat_name);
    sprintf('Computing RS1 structural properties for: %s', mat_name)
    
    % 1. MinNonzerosPerRow
    nnz_per_row = sum(A ~= 0, 2);
    min_non_zeros_per_row = full(min(nnz_per_row));
    sprintf('Minimum no. of non zeros per row: %d', min_non_zeros_per_row)

    % 2. NonZeroPatternSymmetryV1 (Checks the nonzero pattern symmetry, if symmetric returns 1) 
    % if A and A(Transpose) have same nonzero pattern then returns 1)
    non_zero_pattern_A = A ~= 0;
    non_zero_pattern_A_tran = A' ~= 0;
    non_zero_pattern_symmetryV1 = isequal(non_zero_pattern_A, non_zero_pattern_A_tran);
    sprintf('Non zero pattern symmetry: %d', non_zero_pattern_symmetryV1)

    % 3. InfinityNorm (maximum of the absolute row sums)
    inf_norm = norm(A,'inf'); 
    sprintf('Infinity Norm: %d', full(inf_norm))

    % 4. ColumnVariance
    [out_rows, out_cols] = size(A);
    col_variance = 0;
    for i = 1: out_cols
    	col_variance = col_variance + var(A(:,i));
    end
    col_variance = col_variance / out_cols;
    sprintf('Average Column variance: %d\n ', col_variance)
    
	% Another way to do it: just variance of the nonzeros. faster.
	col_variance_nonzeros = sum(arrayfun(@(n) var(nonzeros(A(:,n))), 1:size(A,2)));    
    col_variance_nonzeros = col_variance_nonzeros / out_cols;
    sprintf('Average Column variance of the non zeros: %d\n', col_variance_nonzeros)
     
    % 5. lowerBandwidth
    lower_bw = bandwidth(A, 'lower');
    sprintf('Lower bandwidth %d', lower_bw) 

    % 6. DiagonalNonZeros
    % counts the number of nonzeros on the diagonal
    D = diag(A);
    nnz_diag = nnz(D);
    sprintf('Diagonal non zeros: %d', nnz_diag)
    
    % 7. DiagonalAverage
    % Computes the average of the absolute values of the diagonal elements of a matrix
    diag_avg = mean(abs(diag(A)));
    sprintf('Diagonal average: %d', full(diag_avg))
    
    fprintf(out_file, '%s, %f, %f, %f, %f, %f, %f, %f, %f \n', mat_name, min_non_zeros_per_row, non_zero_pattern_symmetryV1, full(inf_norm), col_variance, col_variance_nonzeros, lower_bw, nnz_diag, full(diag_avg) );
end
