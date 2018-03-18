% Computes the smallest eigen value with the eig() and the power method applied twice.
% Output: Smallest eigen value
% Date: January 29, 2018

%Pseudocode
% This approach uses Power Method twice
%1. Get A
%2. Get lamda_max of A by Power method
%3. B = A - (lamda_max * I) %shift
%4. Get lamda_max of B by Power method
%5. lamda_min of A = (lamda_max of B) + (lamda_max of A) 

file_location = '/Users/kanikas/Documents/research/data_(mat)_matlab_format/mat_files';
all_files = dir(file_location);
all_names = { all_files.name };

out_file = fopen('/Users/kanikas/Documents/MatrixFree/Matrix-free/ConditionNo/eig_val_min_max_iter4_rel_err.csv','a') ;
fprintf(out_file,'Filename, Dimension , No. of Iterations, True Eigen Value, Eigen Value from Power method, Absolute Error, Relative Error, Condition No. from power method, True Condition number, Relative error for condition no. \n');


for i = 4: length(all_names)
    if all_files(i).bytes > 1000000 % less than 1 MB for now
        continue
    else 
        in_file = strcat(file_location,'/', all_files(i).name);
        A = load(in_file);
        sprintf('Computing values for file: %s', all_files(i).name)
        %Step 1: Get A
        %A = [1 1 -2; -1 2 1; 0 1 -1];
        %A = [-1,8,7; 4,3,10; 11,-9,17];
        %A = rand(1000);
        %A = load('/Users/kanikas/Documents/research/data_(mat)_matlab_format/mat_files/2cubes_sphere.mat');
        %A = [ 4 1 0; 0 2 1; 0 0 -1];
        %A = [0 1 0; 0 -2 1; 0 0 -5];
        %A = [4 3 8; 3 2 5; 8 5 6];
        %A = [ -5 -1; 4 -2]; 
        %A = [0 2 14; 2 -8 6; 14 6 9];
        %A = [12 -1; -1 12];
        A = A.Problem.A;
        n = size(A,1);
        w = rand(n,1); % choosing a random vector for initial guess
        w = w/norm(w);
        cntr = 0;

        %Step 2: Get largest eigen value of A by power method
        lamda_max_A = power_method(A,w); 
        I = speye(n); % for very large Identity matrix eye(n) doesnt work: Has error: Error using eye Requested 60000x60000 (26.8GB) array exceeds maximum array size preference. Creation of arrays greater than this limit may take a long time and cause MATLAB to become unresponsive. See array
        %size limit or preference panel for more information.

        %Step 3: Get B = A - (lamda_max * I) %shift
        B = A - (lamda_max_A * I); 

        %Step 4: Get lamda_max of B by Power method
        lamda_max_B = power_method(B,w);

        %Step 5: Get smallest Eigen value lamda_min of A = (lamda_max of B) + (lamda_max of A) 
        lamda_min_A = lamda_max_B + lamda_max_A;

        true_smallest_eigen_val = eigs(A, 1, 'sr');
        sprintf('True smallest eigen value of A: %f ', true_smallest_eigen_val)
        sprintf('Smallest Eigen value from double power method approach: %f ', lamda_min_A)
        abs_err_small = abs(true_smallest_eigen_val) - abs(lamda_min_A);
        abs_err_small = abs(abs_err_small)
        rel_err_small = abs(abs_err_small / true_smallest_eigen_val)
        condn_num = abs(lamda_max_A)/abs(lamda_min_A);
        sprintf('Condition number from double power method approach: %f', condn_num)
        true_condn_num = cond(A);
        sprintf('True Condition number: %f', true_condn_num)
        rel_err_cond = abs(condn_num/true_cond_num);
    end
end

%fprintf(out_file, '%s, %i, %i, %f, %f, %f, %f \n', all_files(i).name, n, iter - 1, true_eig_val, lamda, abs_err, rel_err, condn_num, true_condn_num, rel_err_cond);
fclose(out_file);

function lamda = power_method(A,w) 
    iter = 1;
    %abs_err = 0;
    %rel_err = 0;
    
    while (iter < 5) %fixed # iterations
        %cntr = cntr + 1;
        iter = iter + 1;
        b = A * w; %  matrix-vector product 
        lamda = max(abs(b));
        lamda_with_sign = max(b);
        if lamda == lamda_with_sign
            lamda_sign = 1; % positive
        else   
            lamda_sign = -1; % negative      
        end
        lamda = lamda * lamda_sign;
        w = b/lamda; %x1 = x0/?
    end

    val = max(abs(eigs(A)));
    val_with_sign = max(eigs(A));
    if val == val_with_sign
       val_sign = 1; % positive
    else   
       val_sign = -1; % negative      
    end
       val = val * val_sign;

    %sprintf('Dimension: %i', n)
    true_eig_val = val;
    sprintf('Eigen value: %f', true_eig_val) % largest eigen value
    sprintf('Eigen value from power method: %f', lamda)
    abs_err = abs(true_eig_val) - abs(lamda);
    abs_err = abs(abs_err)
    rel_err = abs(abs_err / true_eig_val)
    %fprintf(out_file, '%s, %i, %i, %i, %i, %i \n', all_files(i).name, n, iter, true_eig_val, lamda, abs(true_eig_val - lamda));
end


