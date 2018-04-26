% Computes various structural properties of a square matrix.
% Output: nnz, max. no. of non zeros per row etc. 
% Date: February 13, 2018

%A = load('/Users/kanikas/Documents/research/data_(mat)_matlab_format/mat_files/494_bus.mat');
%A = load('/Users/kanikas/Documents/DONE/d_ss.mat');

% cd '/Users/kanikas/Documents/DONE'
% file_location = '/Users/kanikas/Documents/DONE';
cd '/Volumes/Kank/MatrixFree/DONE';
file_location = '/Volumes/Kank/MatrixFree/DONE';

all_files = dir(file_location);
all_names = { all_files.name };
global in_file;

out_file = fopen('/Volumes/Kank/MatrixFree/struct_properties_all_matr_final_v4.csv','a') ;
fprintf(out_file,'Matrix name, Dimension, nnz, MaxNonzerosPerRow, MinNonzerosPerRow, AvgNonzerosPerRow, AbsoluteNonZeroSum, Symmetricity, AvgDiagDist, NonZeroPatternSymmetryV1, Trace, AbsoluteTrace, OneNorm, InfinityNorm, FrobeniusNorm, SymmetricInfinityNorm, SymmetricFrobeniusNorm, AntiSymmetricInfinityNorm, AntiSymmetricFrobeniusNorm, DiagonalAverage, RowDiagonalDominance, ColDiagonalDominance, RowVariance, RowVarianceNonZeros, ColumnVariance, ColVarianceNonZeros, lowerBandwidth, upperBandwidth, DiagonalSign, DiagonalNonZeros \n');


for j = 4: length(all_names)
    sprintf('Computing structural properties for: %s', all_files(j).name)
    in_file = strcat(file_location,'/', all_files(j).name);
    A= load(in_file);
    A = A.Problem.A;

    %1. Dimension 
    n = size(A,1);
    sprintf('Matrix dimension: %d', n)

    %2. No. of non-zeros (nnz)
    non_zeros = nnz(A);
    sprintf('No. of non-zeros: %d', non_zeros);

    nnz_per_row = sum(A ~= 0, 2);
    %3. MaxNonzerosPerRow
    max_non_zeros_per_row = full(max(nnz_per_row)); % to convert sparse to dense
    sprintf('Maximum no. of non zeros per row: %d', max_non_zeros_per_row);

    %4. MinNonzerosPerRow
    min_non_zeros_per_row = full(min(nnz_per_row));
    sprintf('Minimum no. of non zeros per row: %d', min_non_zeros_per_row);

    % 5. AvgNonzerosPerRow
    avg_non_zeros_per_row = sum(nnz_per_row)/n;
    sprintf('Average no. of non zeros per row: %d', full(avg_non_zeros_per_row));

    % 6. AbsoluteNonZeroSum (Sum of the absolute values of all the non zeros in matrix) 
    sum_abs_non_zero = sum(sum(abs(A)));
    sprintf('Absolute sum of all non zeros: %d', full(sum_abs_non_zero));

    % 7. NumericValueSymmetryV1/Symmetricity (Checks the numerical symmetry: 1 means symmetric, 0 means non-symmetric)
    symmetricity = issymmetric(A);
    sprintf('Symmetricity: 1 means symmetric, 0 means non-symmetric: %d', symmetricity);

    % 8. NonZeroPatternSymmetryV1 (Checks the nonzero pattern symmetry, if symmetric returns 1) 
    % if A and A(Transpose) have same nonzero pattern then returns 1)
    non_zero_pattern_A = A ~= 0;
    non_zero_pattern_A_tran = A' ~= 0;
    non_zero_pattern_symmetryV1 = isequal(non_zero_pattern_A, non_zero_pattern_A_tran);
    sprintf('Non zero pattern symmetry: %d', non_zero_pattern_symmetryV1);

    % 9. Trace
    tr = sum(diag(A)); %or trace(A)
    sprintf('Trace: %d', full(tr));

    % 10. AbsoluteTrace
    abs_trace = sum(diag(abs(A))); 
    sprintf('Absolute Trace: %d', full(abs_trace));

    % 11. OneNorm (maximum column sum)
    one_norm = max(sum(abs(A)));
    sprintf('One Norm: %d', full(one_norm));

    % 12. InfinityNorm (maximum of the absolute row sums)
    inf_norm = norm(A,'inf'); %max(sum(abs(A,2)));
    sprintf('Infinity Norm: %d', full(inf_norm)); %or norm(A,'inf');

    % 13. FrobeniusNorm
    % square root of the sum of the absolute squares of its elements
    % fro_norm = sqrt(sum(abs(A*A')));
    fro_norm = norm(A, 'fro');
    sprintf('Frobenius Norm: %d', fro_norm) ;

    % 14. SymmetricInfinityNorm
    % Computes the infinity norm of the symmetric part of the matrix
    symm_A = 1/2 * (A + A');
    sym_inf_norm = max(sum(abs(symm_A),2));
    sprintf('Symmetric Infinity Norm: %d', full(sym_inf_norm));

    % 15. SymmetricFrobeniusNorm
    % Frobenius norm of the symmetric part of the matrix
    sym_fro_norm = norm(symm_A,'fro');
    sprintf('Symmetric Frobenius Norm: %d', full(sym_fro_norm));

    % 16. AntiSymmetricInfinityNorm
    anti_symm_A = 1/2 * (A - A');
    anti_sym_inf_norm = max(sum(abs(anti_symm_A),2)); 
    sprintf('Anti symmetric Infinity Norm: %d', full(anti_sym_inf_norm));

    % 17. AntiSymmetricFrobeniusNorm
    anti_sym_fro_norm = norm(anti_symm_A,'fro'); 
    sprintf('Anti symmetric Frobenius Norm: %d', anti_sym_fro_norm);

    % 18. DiagonalAverage
    % Computes the average of the absolute values of the diagonal elements of a matrix
    diag_avg = mean(abs(diag(A)));
    sprintf('Diagonal average: %d', full(diag_avg));

    % 19. Average diagonal distance
    % Average distance of nonzero diagonal to the main diagonal
    % first column with the non-zero
    dist = [find(A(1,:)) find(A(:,1))']; % distance of all the non zero diagonals
    avg_diag_dist = mean(dist-1); % subtract the nonzero diag distance from the main diagonal
    sprintf('Average diagonal distance: %d', avg_diag_dist);

    % 20. RowDiagonalDominance 
    % 21. ColDiagonalDominance
    % computes the (row) diagonal dominance of a matrix
    % returns 0 if it's not
    % returns 1 if it's weakly diagonally dominant
    if issymmetric(A)
        if all((2*abs(diag(A))) >= sum(abs(A),2)) == 0
            row_diag_dominance = 1; 
            col_diag_dominance = 1; 
        else
            row_diag_dominance = 0;
            col_diag_dominance = 0 ;

        end
    else 
        row_diag_dominance = 0;
        col_diag_dominance = 0 ;
    end
    sprintf('Row diagonal dominance: %d', row_diag_dominance);
    sprintf('Column diagonal dominance: %d', col_diag_dominance);

    % 22. RowVariance
    % computes the row variance of a matrix
    row_variance = 0;
    [out_rows, out_cols] = size(A);
    for i = 1: out_rows
            row_variance = row_variance + var(A(i,:));
    end
    row_variance = row_variance / out_rows;
    sprintf('Row variance %d', row_variance);
    
    % Another way to do it: just variance of the nonzeros. faster.
    row_variance_nonzeros = sum(arrayfun(@(n) var(nonzeros(A(n,:))), 1:size(A,1)));
    row_variance_nonzeros = row_variance_nonzeros / out_cols;
    sprintf('Average row variance of the non zeros: %d', row_variance_nonzeros);

    % 23. ColumnVariance
    [out_rows, out_cols] = size(A);
    col_variance = 0;
    for i = 1: out_cols
    	col_variance = col_variance + var(A(:,i));
    end
    col_variance = col_variance / out_cols;
    sprintf('Average Column variance: %d', col_variance);
    
	% Another way to do it: just variance of the nonzeros. faster.
	col_variance_nonzeros = sum(arrayfun(@(n) var(nonzeros(A(:,n))), 1:size(A,2)));    
    col_variance_nonzeros = col_variance_nonzeros / out_cols;
    sprintf('Average Column variance of the non zeros: %d', col_variance_nonzeros)

    % 24. lowerBandwidth
    lower_bw = bandwidth(A, 'lower');
    sprintf('Lower bandwidth %d', lower_bw);

    % 25. upperBandwidth
    %[lower,upper] = bandwidth(A);
    upper_bw = bandwidth(A, 'upper');
    sprintf('Upper bandwidth %d', upper_bw);
  
    % 26. DiagonalSign
    % indicates the diagonal sign pattern
    % -2 all negative, -1 some are negative, 0 all zero, 1 some are positive, 2 all positive
    D = diag(A);
    if all(D(:) == 0)== 1 % all elements are zero
        diag_sign = 0;
    else
        if (all(D(:) < 0)) == 1 % all are negative
            diag_sign = -2;
        elseif (any(D(:) < 0)) == 1 % some are nonpositive
            diag_sign = -1;
        elseif (any(D(:) > 0)) == 1 % some are positive
            diag_sign = 1;
        elseif (all(D(:) > 0)) == 1 % all are positive
            diag_sign = 2;
        end
    end
    sprintf('Diagonal sign: %d', diag_sign);

    % 27. DiagonalNonZeros
    % counts the number of nonzeros on the diagonal
    nnz_diag = nnz(D);
    sprintf('Diagonal non zeros: %d', nnz_diag);
    sprintf('name:: %s', all_files(j).name);

    fprintf(out_file, '%s, %i, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f \n', all_files(j).name, n, non_zeros, max_non_zeros_per_row, min_non_zeros_per_row, full(avg_non_zeros_per_row), full(sum_abs_non_zero), symmetricity, avg_diag_dist, non_zero_pattern_symmetryV1, full(tr), full(abs_trace), full(one_norm), full(inf_norm), fro_norm, full(sym_inf_norm), full(sym_fro_norm), full(anti_sym_inf_norm), anti_sym_fro_norm, full(diag_avg), row_diag_dominance, col_diag_dominance, row_variance, row_variance_nonzeros, col_variance, col_variance_nonzeros, lower_bw, upper_bw, diag_sign, nnz_diag );
    movefile(in_file,'../DONE2')
end
