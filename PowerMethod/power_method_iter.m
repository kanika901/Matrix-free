% Computes the larget eigen value with the eig() and the power method.
% Output: Largest eigen value
% Date: November 18, 2017

A = [1 1 -2; -1 2 1; 0 1 -1];
A = [-1,8,7; 4,3,10; 11,-9,17];
A = rand(1000);

n = size(A,1);
w = rand(n,1); % choosing a random vector for initial guess
w = w/norm(w);
cntr = 0;
iter = 0;
while (iter < 5) %fixed # iterations
    cntr = cntr + 1;
    iter = iter + 1;
	b = A * w; % only matrix-vector product with A % 1
	lamda = max(abs(b));
	w = b/lamda; %x1 = x0/lamda
end
val = eig(A);

sprintf('Dimension: %i', n)
sprintf('Eigen value: %f', val(1)) % largest eigen value
sprintf('Eigen value from power method: %f', lamda)
sprintf('No. of matrix-vector products: %i', cntr) 