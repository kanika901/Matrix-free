% Computes the smallest eigen value with the eig() and the inverse power method.
% Output: Smallest eigen value
% Date: November 23, 2017

A = [-1,8,7; 4,3,10; 11,-9,17];
A = [1 1 -2;-1 2 1; 0 1 -1];
A = [-4 1 1;0 3 1;-2 0 15];
A = [2 2 -1;-5 9 -3; -4 4 1];
A = [3 5 2; 5 12 7; 2 7 5];
A = rand(1000);

n = size(A,1);
y = rand(n,1); % choosing a random vector for initial guess
w = w/norm(w);
cntr = 0;
iter = 0;
b = y;
while (iter < 5) %fixed # iterations
    w = A\b;
    norm_w = norm(w);
    b = w/norm_w;
    lamda = 1/norm_w;
    iter = iter + 1;
end
val = eigs(A,1,'sm');

sprintf('Dimension: %i', n)
sprintf('Smallest Eigen value: %f', val(1)) % smallest eigen value
sprintf('Smallest Eigen value from inverse power method: %f', lamda)
