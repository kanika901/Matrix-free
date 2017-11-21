% Computes the larget eigen value with the eig() and the power method.
% Output: Largest eigen value
% Date: November 8, 2017

A = [1 1 -2;-1 2 1; 0 1 -1];
A = [-1,8,7; 4,3,10; 11,-9,17];
%A = rand(1000);

n = size(A,1);
w = rand(n,1); % choosing a random vector for initial guess
w = w/norm(w);
res = 1;
tol = 1e-8;
cntr = 0;
while (res > tol)
    cntr = cntr + 1;
	b = A*w; % only matrix-vector product with A % 1
	lamda = max(abs(b));
	w = b/lamda; %x1 = x0/lamda
	res = norm(A*w-lamda*w); %x1 -x0 % 2 and 3
end
val = eig(A);
sprintf('Dimension: %i\n', n )
sprintf('Eigen value: %f\n', val(1)) % largest eigen value
sprintf('Eigen value from power method: %f\n', lamda)
