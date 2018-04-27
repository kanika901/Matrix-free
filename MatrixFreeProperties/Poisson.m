% returns the block tridiagonal stencil function of order n^2 resulting
% from discretizing Poisson's equation with the 5-point operator on an
% n-by-n mesh.
function px = Poisson(x, n)
	assert(length(x) == n^2, 'length of x must be n^2');
	px = 4 * x;
	% 1-diagonals
	px(2:end)   = px(2:end) - x(1:end-1);                      % upper
	px(1:end-1) = px(1:end-1) - x(2:end);                      % lower
	% n-diagonals
	px(1:end-n) = px(1:end-n) - x(n+1:end);                    % upper
	px(n+1:end) = px(n+1:end) - x(1:end-n);                    % lower
	% Add back in some 1-diagonals; the poisson matrix is 0 at
	% (k*n, k*n+1) and (k*n+1, k*n) for 1 <= k <= n^2 - n
	px(n:n:end-n) = px(n:n:end-n) + x(n+1:n:end-n+1);          % upper
	px(n+1:n:end-n+1) = px(n+1:n:end-n+1) + x(n:n:end-n);      % lower
end
