% Calculate the jacobian about x at a standard basis vector v.
function Jv = Jv_approx_basis(F, x, v)
	assert(sum(nnz(v)) == 1, ...
			'v must be a standard basis vector; for others use Jv_approx');
	assert(norm(v,1) == 1, 'v must be a standard basis vector');
	% Calculate perturbation. Contrast this with JV_APPROX;
	% We only multiply by the standard basis vectors so we can simplify
	dim = length(x);
	Jv = zeros(dim, 1);
	per = sqrt(eps) * 2 / dim;
	R = F(x); % unperturbed residual
	xper = x + per * v;    % perturbed vector
	Rper = F(xper);        % perturbed residual
	Jv = (Rper - R) / per; % approximation of jacobian on krylov vector
end

