% Author: Sam Pollard
% Code taken from Kanika Sood (structural_properties.m), Cheuk Lau
% (JV_APPROX.m) and Arvind K. Saibaba (subspaceiteration.m),
% and the authors of PETSc
% http://www.mcs.anl.gov/petsc/petsc-current/src/snes/examples/tutorials/ex5s.c.html
path(path,'petsc');

% Set up parameters for the PDEs
tol = 1e-6; % For comparing Jacobians

% Bratu
mx  = 5; % Grid points in the x-direction
my  = 5; % Grid points in the y-direction
lambda = 2.5; % Nonlinearity parameter, 0 <= lambda <= 6.81
run('petsc/ex5_x0.m');
x0B = Vec_ex5_0;
B = @(x) Bratu(x, lambda, mx, my);

% Poisson
p_n = 5; % Poisson (we only support a square grid)
x0P = rand(p_n * p_n, 1);
P = @(x) Poisson(x, p_n);
assert(p_n^2 == mx*my, 'Poisson and Bratu must be the same problem size');

% DrivenCavity
mxD = 5;
myD = 5;
lidvelocity = 10;
grashof = 10^5;
prandtl = 1;
run('petsc/ex19_x0.m');
x0 = Vec_ex19_0;
x0u     = x0(1:4:mxD*myD*4-3);
x0v     = x0(2:4:mxD*myD*4-2);
x0omega = x0(3:4:mxD*myD*4-1);
x0temp  = x0(4:4:mxD*myD*4);
x0D = [x0u, x0v, x0omega, x0temp];
clear x0

% Sanity check: Actually form the Jacobians about x
J_B = FormJacobian(B, x0B);
J_P = FormJacobian(P, x0P);
%A   = gallery('normaldata', n*n, n); % Can also test with random matrix
% This is how you reshape for DrivenCavity as I wrote it.
% This is probably not the best way to do it, since in PETSc if there are
% multiple values, they all get crammed into one Jacobian. This makes sense
% because one property (e.g. u velocity) may affect another (e.g. temp)
D = @(xu, xv, xomega, xtemp) ...
	DrivenCavity(xu, xv, xomega, xtemp, mxD, myD, lidvelocity, grashof, prandtl);
J_D = FormMultiJacobian(D, x0D, mxD, myD);
J_Du     = J_D(:,:,1);
J_Dv     = J_D(:,:,2);
J_Domega = J_D(:,:,3);
J_Dtemp  = J_D(:,:,4);

% Verify Jacobians are correct
% We verify via the following:
% Poisson: Just compare J_P with gallery('poisson',n)
% Bratu: Compare reshape(B(x), mx, my)
%        with ex5m(reshape(x,mx,my), lambda)
% TODO: DrivenCavity. Tricky since we don't have a baseline for correctness
assert(norm(J_P - gallery('poisson', p_n)) < tol, ...
	'Poisson Jacobian is incorrect');
bratu_jac_fn = 'petsc/bratu_100x100_lambda5pt5.mat';
if exist(bratu_jac_fn, 'file') == 2
	load(bratu_jac_fn);
	J_B_KSP = bratu_100x100_lambda5p5;
	J_B_M = FormJacobian(@(x) ...
		reshape(ex5m(reshape(x,mx,my), lambda), [], 1), x0B);
	assert(norm(J_B - J_B_M) < tol, 'Bratu Jacobian is incorrect');
else
	fprintf('Could not find %s\n', bratu_jac_fn);
end
fprintf('Warning: DrivenCavity.m is not verified to be correct\n');

% Calculate matrix properties
[full_propsB, full_nameB] = matrix_full_properties(J_B);
[full_propsP, full_nameP] = matrix_full_properties(J_P);

[free_propsB, free_nameB] = matrix_free_properties(B, x0B);
[free_propsP, free_nameP] = matrix_free_properties(P, x0P);

computed = isfinite(free_propsP);
names = free_nameP(computed)
fullthenfree = [full_propsP(computed); free_propsP(computed)]

% Helper functions
function J = FormJacobian(F, x)
	% Calculate perturbation. Contrast this with JV_APPROX;
	% We only multiply by the standard basis vectors so we can simplify
	dim = length(x);
	J = zeros(dim);
	per = sqrt(eps) * 2 / dim;
	R = F(x); % unperturbed residual
	for i = 1:dim
		v = sparse(i, 1, 1, dim, 1, 1);
		% Approximate J*h_i
		xper = x + per * v; % perturbed vector
		Rper = F(xper); % perturbed residual
		y = (Rper - R) / per; % approximation of jacobian on krylov vector
		J(:,i) = y;
	end
end

function J = FormMultiJacobian(F, x, mx, my)
	% Form the Jacobian with respect to multiple values
	% x is a matrix with the columns corresponding to each input vector
	% Returns J which is a 3-tensor where J(:,:,i) corresponds to x(:,i)
	dim = size(x);
	% dim(1) = mx*my, dim(2) = number of arguments of F
	J = zeros(dim(1),dim(1),dim(2));
	per = sqrt(eps) * 2 / dim(1);
	xsq = zeros(mx,my,dim(2));
	for i = 1:dim(2)
		xsq(:,:,i) = reshape(x(:,i), mx, my);
	end
	R = F(xsq(:,:,1), xsq(:,:,2), xsq(:,:,3), xsq(:,:,4)); % unperturbed residual
	for j = 1:dim(1)
		v = sparse(j, 1, 1, dim(1), 1, 1);
		xper = x + repmat(per * v, 1, dim(2)); % perturbed vector
		Rper = F(reshape(xper(:,1),mx,my), ...
			     reshape(xper(:,2),mx,my), ...
				 reshape(xper(:,3),mx,my), ...
				 reshape(xper(:,4),mx,my)); % perturbed residual
		for i = 1:dim(2)
			R_i = R(:,i);
			Rper_i = Rper(:,i);
			% Approximate J*h_i
			y = (Rper_i - R_i) / per; % approximation of jacobian on krylov vector
			J(:,j,i) = y;
		end
	end
end

