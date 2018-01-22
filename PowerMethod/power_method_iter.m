% Computes the larget eigen value with the eig() and the power method.
% Output: Largest eigen value
% Date: November 18, 2017
%Todo: give dominant eigen value instead of largest

A = [1 1 -2; -1 2 1; 0 1 -1];
A = [-1,8,7; 4,3,10; 11,-9,17];
A = rand(1000);
A = load('/Users/kanikas/Documents/research/data_(mat)_matlab_format/mat_files/2cubes_sphere.mat');
A = [ 4 1 0; 0 2 1; 0 0 -1];
%A = [0 1 0; 0 -2 1; 0 0 -5];
%A = [4 3 8; 3 2 5; 8 5 6];
A = [ -5 -1; 4 -2]; % not working right now because this program should give the dominant eigen value, it gives largest eigen value as of now
%A = [0 2 14; 2 -8 6; 14 6 9];
%A = [12 -1; -1 12];

n = size(A,1);
w = rand(n,1); % choosing a random vector for initial guess
w = w/norm(w);
fprintf('w:: \n')
disp(w)
cntr = 0;
iter = 0;
while (iter < 4) %fixed # iterations
    cntr = cntr + 1;
    iter = iter + 1;
    fprintf('b:: \n')
	b = A * w; %  matrix-vector product 
	lamda = max(abs(b));
    disp(b)
    disp(abs(b))
    disp(max(abs(b)))
    lamda_with_sign = max(b);
    if lamda == lamda_with_sign
        lamda_sign = 1; % positive
    else   
        lamda_sign = -1; % negative      
    end
    lamda = lamda * lamda_sign;
	w = b/lamda; %x1 = x0/?
end

val = max(abs(eig(A)));
val_with_sign = max(eig(A));
if val == val_with_sign
        val_sign = 1; % positive
    else   
        val_sign = -1; % negative      
end
val = val * val_sign;

sprintf('Dimension: %i', n)
sprintf('Eigen value: %f', val) % largest eigen value
sprintf('Eigen value from power method: %f', lamda)
sprintf('No. of matrix-vector products: %i', cntr) 