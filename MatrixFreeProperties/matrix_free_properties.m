% Given a vector function F and a vector x, computes or approximates the
% structural properties of the Jacobian at x only through approximating
% matrix-vector products J_x * v.

% Throughout, J is the jacobian at the vector x and 
% e_i are the standard basis vectors; e.g. e_3 is [0; 0; 1; 0; ... ].

% These approximations are mostly done via computing O(1) columns of J
% then using these as a sample. If rows are needed, we take just horizontal
% slices of the columns (thus an even smaller sample).
% Note that if we know something about F, many of these properties are
% easy (e.g. we know the stencil may be diagonally dominant, or
% nonsymmetric a priori).
% I rank these from 0 to 10, where 10 is it's definitely correct and 0
% means I have no faith in the approximation.

% Some things to keep in mind:
% 0. J*e_i(j) gives J_{ji} and J*e_j(i) gives J_{ij}
% 1. If the system is banded, then we could
%    encode multiple columns of J inside a single matrix-vector product.
%    For example, if we know the maximum diagonal distance is 5, then
%    we could potentially get n/10 columns per J*v by choosing
%    v = e_1 + e_11 + e_21 + ... + e_n
% 2. If the system is symmetric we have things a lot easier since we can
%    compute both a row and a column of a matrix


function [props, names] = matrix_free_properties(F, x)
	names = {'Dimension', 'nnz', 'MaxNonzerosPerRow', 'MinNonzerosPerRow', 'AvgNonzerosPerRow', 'AbsoluteNonZeroSum', 'symmetricity', 'NonZeroPatternSymmetryV1', 'Trace', 'AbsoluteTrace', 'OneNorm', 'InfinityNorm', 'FrobeniusNorm', 'SymmetricInfinityNorm', 'SymmetricFrobeniusNorm', 'AntiSymmetricInfinityNorm', 'AntiSymmetricFrobeniusNorm', 'DiagonalAverage', 'AvgDiagDistance', 'RowDiagonalDominance', 'ColDiagonalDominance', 'RowVariance', 'ColumnVariance', 'lowerBandwidth', 'upperBandwidth', 'DiagonalMean', 'DiagonalSign', 'DiagonalNonZeros'};
	% Parameters that simplify calculations
	% NOTE: Later, we should be able to take in parameters from the user
	% if s/he knows something about F. Symmetry is the most significant.
	symmetric = false; % This is alwys false for now. See note for NumericValueSymmetryV1 

	% Others:
	% # of points on the stencil

	%1. Dimension - 10/10
	n = length(x);
	fprintf('Matrix dimension: %d\n', n);

	% Pick values that probably won't be boundary positions
	% We use these throughout the code
	% Depending on the dimensionality of the data, we may want to
	% choose other points offset of mid
	mid = ceil(n/2);
	Jcols_mid = [mid-2; mid-1; mid; mid+1; mid+2; mid+3; ...
			mid+1-floor(n/5); floor(n/2)+floor(n/5); ...
			mid+1-floor(n/7); floor(n/2)+floor(n/7)];
	Jcols_mid = unique(Jcols_mid);
	n_samples_mid = length(Jcols_mid);
	vmid  = sparse(n, n_samples_mid);
	Jvmid = sparse(n, n_samples_mid);
	for c = 1:length(Jcols_mid)
		vmid(:,c) = sparse(Jcols_mid(c), 1, 1, n, 1, 1);
		Jvmid(:,c) = Jv_approx_basis(F, x, vmid(:,c));
	end

	% Pick values that are probably boundary conditions
	% These are used throughout for guessing minimum values
	Jcols_edge  = [1; 2; 3; 4; 5; n-4; n-3; n-2; n-1; n];
	Jcols_edge  = unique(Jcols_edge);
	n_samples_edge  = length(Jcols_edge);
	vedges  = sparse(n, n_samples_edge);
	Jvedges = sparse(n, n_samples_edge);
	for c = 1:length(Jcols_edge)
		vedges(:,c) = sparse(Jcols_edge(c), 1, 1, n, 1, 1);
		Jvedges(:,c) = Jv_approx_basis(F, x, vedges(:,c));
	end

	% Combine boundary columns with edge columns
	[uniq_e,uniq_idx] = setdiff(Jcols_edge, Jcols_mid);
	Jcols = [Jcols_mid; uniq_e];
	n_samples = length(Jcols);
	Jvsamples = [Jvmid, Jvedges(:,uniq_idx)];
	[Jcols, SI] = sort(Jcols);
	Jvsamples = Jvsamples(:,SI);


	%2. No. of non-zeros (nnz) - 8/10
	sample = nnz(Jvmid);
	non_zeros = n * sample / n_samples_mid;
	fprintf('No. of non-zeros: %d\n', non_zeros);

	%3. MaxNonzerosPerRow - 8/10
	if symmetric
		max_non_zeros_per_row = full(max(sum(Jvmid~=0)));
	else
		rows = Jvmid(Jcols_mid,:);
		max_non_zeros_per_row = full(max(sum(rows~=0)));
	end
	fprintf('Maximum no. of non zeros per row: %f\n', max_non_zeros_per_row)

	%4. MinNonzerosPerRow - 8/10
	if symmetric
		min_non_zeros_per_row = full(min(sum(Jvedges~=0)));
	else
		rows = Jvedges(Jcols_edge,:);
		min_non_zeros_per_row = full(min(sum(rows~=0)));
	end
	fprintf('Minimum no. of non zeros per row: %f\n', min_non_zeros_per_row);

	% 5. AvgNonzerosPerRow - 8/10
	avg_non_zeros_per_row = full(non_zeros) / n;
	fprintf('Average no. of non zeros per row: %f\n', avg_non_zeros_per_row);

	% 6. AbsoluteNonZeroSum (Sum of the absolute values of all the non zeros in matrix) - 8/10
	sum_abs_non_zero = full(sum(sum(abs(Jvmid))))/n_samples_mid * n;
	fprintf('Absolute sum of all non zeros: %d\n', sum_abs_non_zero);

	% 7. NumericValueSymmetryV1 (Checks the numerical symmetry: 1 means symmetric, 0 means non-symmetric) - 3/10
	% symmetry
	% XXX: This is dangerous, since some iterative methods require symmetry.
	%      Eventually, this should be made into a flag or something that gets
	%      passed into the function. Ideally the user of the matrix-free method
	%      nnow if F produces a symmetric jacobian and can thus set 'symmetric'
	%      to true. For now, we do our best guess.
	if symmetric
		symmetricity = 1;
	else
		% For each column we've sampled, check the corresponding row
		symmetricity = 1;
		for j = Jcols'
			symmetricity = symmetricity && ...
					(all(full(Jvsamples(Jcols,j)) == full(Jvsamples(j,:)')));
		end
	end
	fprintf('Symmetricity: 1 means symmetric, 0 means non-symmetric: %d\n', symmetricity);

	% 8. NonZeroPatternSymmetryV1 (Checks the nonzero pattern symmetry, if symmetric returns 1)
	if symmetric || symmetricity
		non_zero_pattern_symmetryV1  = 1;
	else
		non_zero_pattern_symmetryV1  = 1;
		for j = Jcols'
			non_zero_pattern_symmetryV1  = non_zero_pattern_symmetryV1 && ...
					(all(logical(Jvsamples(Jcols,j)) == logical(Jvsamples(j,:)')));
		end
	end
	fprintf('Non zero pattern symmetry: %d\n', non_zero_pattern_symmetryV1);

	% 9. Trace
	% Try 1: tends to underestimate it for these typical stencil problems - 4/10
	% This is from "Randomized Matrix-free Trace and Log-Determinant Estimators"
	% Arvind K. Saibaba, Alen Alexanderian, Ilse C.F. Ipsen
	if symmetric
		p = min(10, ceil(log2(n)));
		q = 1;
		Omega_g = randn(n,p); % Random Gaussian starting guess
		% TODO: See if we can improve Omega_g starting guesses
		tr = randtrace(F,Omega_g,q);
	% Try 2: I think this is biased, but tends to be more accurate - 6/10
	else
		tr = 0;
		for c = 1:length(Jcols)
			r = Jcols(c);
			tr = tr + Jvsamples(r,c);
		end
		tr = n * tr / n_samples;
	end
	fprintf('Trace: %f\n', tr);

	% 10. AbsoluteTrace - 6/10
	abs_trace = 0;
	for c = 1:length(Jcols)
		r = Jcols(c);
		abs_trace = abs_trace + abs(Jvsamples(r,c));
	end
	abs_trace = n * abs_trace / n_samples;
	fprintf('Absolute Trace: %d\n', abs_trace);

	% 11. OneNorm (maximum column sum) - 7/10
	one_norm = full(max(sum(abs(Jvsamples))));
	fprintf('One Norm: %d\n', one_norm);

	% 12. InfinityNorm (maximum of the absolute row sums) - 7/10
	if symmetric
		inf_norm = one_norm;
	else
		rows = Jvsamples(Jcols,:);
		inf_norm = full(max(sum(abs(rows), 2)));
	end
	fprintf('Infinity Norm: %d\n', inf_norm);

	% 13. FrobeniusNorm - 5/10
	% square root of the sum of the absolute squares of its elements
	% This tends to overestimate since the center columns are typically denser
	sample = full(sqrt(sum(sum(abs(Jvmid).^2))));
	fro_norm = n * sample / n_samples_mid;
	fprintf('Frobenius Norm: %d\n', fro_norm);

	% 14. SymmetricInfinityNorm - 6/10
	% Computes the infinity norm of the symmetric part of the matrix
	%symm_A = 1/2 * (A + A');
	if symmetric
		sym_inf_norm = inf_norm;
	else
		rows = Jvsamples(Jcols,:);
		% rows = Jvmid(Jcols_mid,:); % Useful for testing since for n=4 it doesn't have all the columns
		sym_sample = 1/2 * (rows + rows');
		sym_inf_norm = full(max(sum(abs(sym_sample), 2)));
	end
	fprintf('Symmetric Infinity Norm: %d\n', sym_inf_norm);

	% 15. SymmetricFrobeniusNorm - 6/10
	% Frobenius norm of the symmetric part of the matrix
	if symmetric
		sym_fro_norm = fro_norm;
	else
		rows = Jvmid(Jcols_mid,:);
		sym_sample = 1/2 * (rows + rows');
		sample = full(sqrt(sum(sum(abs(sym_sample).^2))));
		sym_fro_norm = n * sample / n_samples_mid;
	end
	fprintf('Symmetric Frobenius Norm: %d\n', sym_fro_norm);

	% 16. AntiSymmetricInfinityNorm - 6/10
	%anti_symm_A = 1/2 * (A - A');
	if symmetric
		anti_sym_inf_norm = 0;
	else
		rows = Jvsamples(Jcols,:);
		sample = 1/2 * (rows - rows');
		anti_sym_inf_norm = full(max(sum(abs(sample), 2)));
	end
	fprintf('Anti symmetric Infinity Norm: %d\n', anti_sym_inf_norm);

	% 17. AntiSymmetricFrobeniusNorm - 6/10
	if symmetric
		anti_sym_fro_norm = 0;
	else
		rows = Jvmid(Jcols_mid,:);
		sym_sample = 1/2 * (rows - rows');
		sample = full(sqrt(sum(sum(abs(sym_sample).^2))));
		anti_sym_fro_norm = n * sample / n_samples_mid;
	end
	fprintf('Anti symmetric Frobenius Norm: %d\n', anti_sym_fro_norm);

	% 18. DiagonalAverage 7/10
	% NOTE: If we can assume the diagonal is stored, in PETSc this would be 10/10
	% Computes the average of the absolute values of the diagonal elements of a matrix
	diag_avg = 0;
	for c = 1:length(Jcols_mid)
		r = Jcols_mid(c);
		diag_avg = diag_avg + abs(Jvmid(r,c));
	end
	diag_avg = diag_avg / n_samples_mid;
	fprintf('Diagonal average: %d\n', diag_avg);

	% 19. AvgDiagDistance - 7/10
	% Average distance of nonzeros to the diagonal
	avg_diag_distance = 0;
	denom = 0;
	for c = 1:size(Jvmid,2)
		nnzlocs = find(Jvmid(:,c));
		if length(nnzlocs) == 0
			continue
		else
			denom = denom + length(nnzlocs) - 1;
		end
		avg_diag_distance = avg_diag_distance + sum(abs(nnzlocs - Jcols_mid(c)));
	end
	avg_diag_distance = avg_diag_distance / denom;
    fprintf('Average diagonal distance: %d\n', full(avg_diag_distance)) 

	% 21. ColDiagonalDominance - 3/10
	% We compute column diagonal dominance first because it's more
	% straightforward, and in the case of symmetry implies RowDiagDominance
	% XXX: This may be dangerous because certain iterative methods such as
	%      Jacobi or Gauss-Seidel require a strictly diagonally dominant
	%      matrix. This only says for sure whether the matrix _isnt_
	%      diagonally dominant but cannot prove the converse with a
	%      constant number of columns of the Jacobian.
	%      See also note for symmetry
	col_diag_dominance = true;
	for c = 1:length(Jcols)
		colsum = sum(abs(Jvsamples(:,c)));
		r = Jcols(c);
		d = abs(Jvsamples(r,c));
		col_diag_dominance = col_diag_dominance && (colsum - 2*d) <= 0;
	end

	% 20. RowDiagonalDominance 3/10
	% computes the (row) diagonal dominance of a matrix
	% XXX: This is an interesting problem.  "How can you determine if the
	%      Jacobian of a function is diagonally dominant?"
	% returns 0 if it's not
	% returns 1 if it's weakly (or strictly) diagonally dominant
	if symmetric
		row_diag_dominance = col_diag_dominance;
	else
		row_diag_dominance = 1;
		rows = Jvsamples(Jcols,:);
		for r = 1:length(Jcols)
			d = abs(rows(r,r));
			row_diag_dominance =  row_diag_dominance && ...
				(sum(abs(rows(r,:))) - 2*d) <= 0;
		end
	end
	fprintf('Row diagonal dominance: %f\n', row_diag_dominance)
	fprintf('Column diagonal dominance: %f\n', col_diag_dominance)

	% 23. ColumnVariance - 6/10
	% We compute column variance first for the same reasons as above
	sample = var(Jvmid);
	col_variance = full(mean(sample));
	fprintf('Column variance %d\n', col_variance);

	% 22. RowVariance - 6/10
	% computes the row variance of a matrix
	if symmetric
		row_variance = col_variance;
	else
		rows = Jvmid(Jcols_mid,:);
		sample = var(rows);
		row_variance = full(mean(sample));
	end
	fprintf('Row variance %d\n', row_variance);

	% In general, I don't think this will work, since the jacobian where
	% F evaluates to 0 is not necessarily 0. However, I think for stencil computations
	% this is true. Also, it relies only on the function F and not a finite-difference
	% approximation of the Jacobian at x, which could be more efficient to explore
	% 24. lowerBandwidth - ?/10
	% 25. upperBandwidth - ?/10
	Fvmid = sparse(n, n_samples_mid);
	for c = 1:n_samples_mid
		Fvmid(:,c) = F(vmid(:,c));
	end
	upper_bw = 0;
	lower_bw = n;
	for c = 1:length(Jcols_mid)
		r = Jcols_mid(c);
		upper_bw = max(upper_bw, r - min(find(Fvmid(:,c))));
		lower_bw = min(lower_bw, max(find(Fvmid(:,c))) - r);
	end

	fprintf('Lower bandwidth %d\n', lower_bw);
	fprintf('Upper bandwidth %d\n', upper_bw);

	% I skip these because I'm not sure how Petsc will end up storing the
	% Jacobian in matrix free. If this is not stored, then you will have
	% to do the trick you find elsewhere with rows = Jvmid(Jcols_mid,:),
	% and you get a limited number of diagonal entries with diag(rows).
	% and then get
	% 26. DiagonalMean
	%D = diag(A);
	%[diag_rows, c] = size(D);
	%diag_mean = 0.0;
	%zero_count = 0;

	%for i = 1:diag_rows
	%	 diag_mean = diag_mean +  D(i);
	%	 if D(i) == 0
	%		 zero_count = zero_count + 1;
	%	 end
	%end
	%diag_mean = diag_mean / diag_rows;
	diag_mean = nan;
	fprintf('Diagonal mean: %d\n', diag_mean);

	% 27. DiagonalSign
	%// indicates the diagonal sign pattern
	%// -2 all negative, -1 nonpositive, 0 all zero, 1 nonnegative, 2 all positive,
	%// 3 some negative,some or no zero,some positive
	%if all(D(:) == 0)== 1 % all elements are zero
	%	 diag_sign = 0;
	%else
	%	 if (all(D(:) < 0)) == 1 % all are negative
	%		 diag_sign = -2;
	%	 elseif (any(D(:) < 0)) == 1 % some are nonpositive
	%		 diag_sign = -1;
	%	 elseif (any(D(:) > 0)) == 1 % some are positive
	%		 diag_sign = -1;
	%	 elseif (all(D(:) > 0)) == 1 % all are positive
	%		 diag_sign = 2;
	%	 end
	%end
	diag_sign = nan;
	fprintf('Diagonal sign: %d\n', diag_sign);

	% 28. DiagonalNonZeros
	% counts the number of nonzeros on the diagonal
	%nnz_diag = nnz(D);
	nnz_diag = nan;
	fprintf('Diagonal non zeros: %d\n', nnz_diag);

	props = [ ...
		n, non_zeros, max_non_zeros_per_row, min_non_zeros_per_row, ...
		full(avg_non_zeros_per_row), full(sum_abs_non_zero), symmetricity,       ...
		non_zero_pattern_symmetryV1, full(tr), full(abs_trace), full(one_norm),  ...
		full(inf_norm), fro_norm, full(sym_inf_norm), full(sym_fro_norm),        ...
		full(anti_sym_inf_norm), anti_sym_fro_norm, full(diag_avg),              ...
		avg_diag_distance, row_diag_dominance, col_diag_dominance, row_variance, ...
		col_variance, lower_bw, upper_bw, diag_mean, diag_sign, nnz_diag];
end
