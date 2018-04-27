% Computes various structural properties of a square matrix.
% Output: nnz, max. no. of non zeros per row etc. 
% Date: March 11, 2018
% Author: Kanika Sood
% Modified by Sam Pollard to make it into a function

function [props, names] = matrix_full_properties(A)
	names = 'Dimension, nnz, MaxNonzerosPerRow, MinNonzerosPerRow, AvgNonzerosPerRow, AbsoluteNonZeroSum, symmetricity, NonZeroPatternSymmetryV1, Trace, AbsoluteTrace, OneNorm, InfinityNorm, FrobeniusNorm, SymmetricInfinityNorm, SymmetricFrobeniusNorm, AntiSymmetricInfinityNorm, AntiSymmetricFrobeniusNorm, DiagonalAverage, AvgDiagDistance, RowDiagonalDominance, ColDiagonalDominance, RowVariance, ColumnVariance, lowerBandwidth, upperBandwidth, DiagonalMean, DiagonalSign, DiagonalNonZeros'; 

    %1. Dimension 
    n = size(A,1);
    fprintf('Matrix dimension: %d\n', n);

    %2. No. of non-zeros (nnz)
	non_zeros = nnz(A);
    fprintf('No. of non-zeros: %d\n', non_zeros);

    nnz_per_row = sum(A ~= 0, 2);
    %3. MaxNonzerosPerRow
    max_non_zeros_per_row = full(max(nnz_per_row)); % to convert sparse to dense
    fprintf('Maximum no. of non zeros per row: %f\n', max_non_zeros_per_row)

    %4. MinNonzerosPerRow
    min_non_zeros_per_row = full(min(nnz_per_row));
    fprintf('Minimum no. of non zeros per row: %f\n', min_non_zeros_per_row)

    % 5. AvgNonzerosPerRow
    avg_non_zeros_per_row = sum(nnz_per_row)/n;
    fprintf('Average no. of non zeros per row: %f\n', full(avg_non_zeros_per_row))

    % 6. AbsoluteNonZeroSum (Sum of the absolute values of all the non zeros in matrix) 
    sum_abs_non_zero = sum(sum(abs(A)));
    fprintf('Absolute sum of all non zeros: %d\n', full(sum_abs_non_zero))

    % 7. NumericValueSymmetryV1 (Checks the numerical symmetry: 1 means symmetric, 0 means non-symmetric)
    symmetricity = issymmetric(A);
    fprintf('Symmetricity: 1 means symmetric, 0 means non-symmetric: %d\n', symmetricity);

    % 8. NonZeroPatternSymmetryV1 (Checks the nonzero pattern symmetry, if symmetric returns 1) 
    % if A and A(Transpose) have same nonzero pattern then returns 1)
    non_zero_pattern_A = A ~= 0;
    non_zero_pattern_A_tran = A' ~= 0;
    non_zero_pattern_symmetryV1 = isequal(non_zero_pattern_A, non_zero_pattern_A_tran);
    fprintf('Non zero pattern symmetry: %d\n', non_zero_pattern_symmetryV1)

    % 9. Trace
    tr = sum(diag(A)); %or trace(A)
    fprintf('Trace: %d\n', full(tr))

    % 10. AbsoluteTrace
    abs_trace = sum(diag(abs(A))); 
    fprintf('Absolute Trace: %d\n', full(abs_trace))

    % 11. OneNorm (maximum column sum)
    one_norm = max(sum(abs(A)));
    fprintf('One Norm: %d\n', full(one_norm))

    % 12. InfinityNorm (maximum of the absolute row sums)
    inf_norm = norm(A,'inf'); %max(abs(sum(A,2)));
    fprintf('Infinity Norm: %d\n', full(inf_norm)) %or norm(A,'inf')

    % 13. FrobeniusNorm
    % square root of the sum of the absolute squares of its elements
    %fro_norm = sqrt(sum(abs(A*A')));
    fro_norm = norm(A, 'fro');
    fprintf('Frobenius Norm: %d\n', fro_norm) 

    % 14. SymmetricInfinityNorm
    % Computes the infinity norm of the symmetric part of the matrix
    symm_A = 1/2 * (A + A');
    sym_inf_norm = max(sum(abs(symm_A),2)); % test it
    fprintf('Symmetric Infinity Norm: %d\n', full(sym_inf_norm))

    % 15. SymmetricFrobeniusNorm
    % Frobenius norm of the symmetric part of the matrix
    sym_fro_norm = norm(symm_A,'fro');
    fprintf('Symmetric Frobenius Norm: %d\n', full(sym_fro_norm)) 

    % 16. AntiSymmetricInfinityNorm
    anti_symm_A = 1/2 * (A - A');
    anti_sym_inf_norm = max(sum(abs(anti_symm_A),2)); % test it
    fprintf('Anti symmetric Infinity Norm: %d\n', full(anti_sym_inf_norm))

    % 17. AntiSymmetricFrobeniusNorm
    anti_sym_fro_norm = norm(anti_symm_A,'fro'); % test it
    fprintf('Anti symmetric Frobenius Norm: %d\n', anti_sym_fro_norm) 

    % 18. DiagonalAverage
    % Computes the average of the absolute values of the diagonal elements of a matrix
    diag_avg = mean(abs(diag(A)));
    fprintf('Diagonal average: %d\n', full(diag_avg)) 


    % 19. AvgDiagDistance
    % Average distance of nonzeros to the diagonal
	avg_diag_distance = 0;
	denom = 0;
	for c = 1:size(A,2)
		nnzlocs = find(A(:,c));
		if length(nnzlocs) == 0
			continue
		else
			denom = denom + length(nnzlocs) - 1;
		end
		avg_diag_distance = avg_diag_distance + sum(abs(nnzlocs - c));
	end
	avg_diag_distance = avg_diag_distance / denom;
    fprintf('Average diagonal distance: %d\n', full(avg_diag_distance)) 

    % 20. RowDiagonalDominance 
    % 21. ColDiagonalDominance
    % computes the (row) diagonal dominance of a matrix
    % returns 0 if it's not
    % returns 1 if it's weakly (or strictly) diagonally dominant

    if issymmetric(A)
        if all((2*abs(diag(A))) >= sum(abs(A),2)) == 0
            row_diag_dominance = 1; 
            col_diag_dominance = 1; 
        else
            row_diag_dominance = 0;
            col_diag_dominance = 0;

        end
    else 
        row_diag_dominance = 0;
        col_diag_dominance = 0;
    end
    fprintf('Row diagonal dominance: %f\n', row_diag_dominance)
    fprintf('Column diagonal dominance: %f\n', col_diag_dominance)

    % 22. RowVariance
    % computes the row variance of a matrix
    %[out_rows, out_cols] = size(A);
	%row_variance = 0;
    %for i = 1: out_rows
	%	row_variance = row_variance + var(A(i,:));
    %end

	% Another way to do it: just variance of the nonzeros. faster.
	row_variance = sum(arrayfun(@(n) var(nonzeros(A(n,:))), 1:size(A,1)));

    row_variance = row_variance / size(A,1);
    fprintf('Row variance %d\n', row_variance)

    % 23. ColumnVariance
	%col_variance = 0;
    %for i = 1: out_cols
    %	col_variance = col_variance + var(A(:,i));
    %end

	% Another way to do it: just variance of the nonzeros. faster.
	col_variance = sum(arrayfun(@(n) var(nonzeros(A(:,n))), 1:size(A,2)));

    col_variance = col_variance / size(A,2);
    fprintf('Column variance %d\n', col_variance)

    % 24. lowerBandwidth
    lower_bw = bandwidth(A, 'lower');
    fprintf('Lower bandwidth %d\n', lower_bw)

    % 25. upperBandwidth
    %[lower,upper] = bandwidth(A);
    upper_bw = bandwidth(A, 'upper');
    fprintf('Upper bandwidth %d\n', upper_bw)

    % 26. DiagonalMean
    D = diag(A);
    [diag_rows, c] = size(D);
    diag_mean = 0.0;
    zero_count = 0;

    for i = 1:diag_rows
        diag_mean = diag_mean +  D(i);
        if D(i) == 0
            zero_count = zero_count + 1;
        end
    end
    diag_mean = diag_mean / diag_rows;
    fprintf('Diagonal mean: %d ', diag_mean)    

    % 27. DiagonalSign
    %// indicates the diagonal sign pattern
    %// -2 all negative, -1 nonpositive, 0 all zero, 1 nonnegative, 2 all positive, 
    %// 3 some negative,some or no zero,some positive
    if all(D(:) == 0)== 1 % all elements are zero
        diag_sign = 0;
    else
        if (all(D(:) < 0)) == 1 % all are negative
            diag_sign = -2;
        elseif (any(D(:) < 0)) == 1 % some are nonpositive
            diag_sign = -1;
        elseif (any(D(:) > 0)) == 1 % some are positive
            diag_sign = -1;
        elseif (all(D(:) > 0)) == 1 % all are positive
            diag_sign = 2;
        end
    end
    fprintf('Diagonal sign: %d\n', diag_sign)

    % 28. DiagonalNonZeros
    % counts the number of nonzeros on the diagonal
    nnz_diag = nnz(D);
    fprintf('Diagonal non zeros: %d\n', nnz_diag)

	props = [ ...
		n, non_zeros, max_non_zeros_per_row, min_non_zeros_per_row, ...
		full(avg_non_zeros_per_row), full(sum_abs_non_zero), symmetricity,       ...
		non_zero_pattern_symmetryV1, full(tr), full(abs_trace), full(one_norm),  ...
		full(inf_norm), fro_norm, full(sym_inf_norm), full(sym_fro_norm),        ...
		full(anti_sym_inf_norm), anti_sym_fro_norm, full(diag_avg),              ...
		avg_diag_distance, row_diag_dominance, col_diag_dominance, row_variance, ...
		col_variance, lower_bw, upper_bw, diag_mean, diag_sign, nnz_diag];
end

