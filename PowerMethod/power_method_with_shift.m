% Computes the smallest eigen value with the eig() and the power method applied twice.
% Output: Smallest eigen value
% Date: January 26, 2018

%Pseudocode
% This approach uses Power Method twice
%1. Get A
%2. Get lamda_max of A by Power method
%3. B = A - (lamda_max * I) %shift
%4. Get lamda_max of B by Power method
%5. lamda_min of A = (lamda_max of B) + (lamda_max of A) 