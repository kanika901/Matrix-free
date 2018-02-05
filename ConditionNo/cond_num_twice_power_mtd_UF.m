% Computes the smallest eigen value with the eig() and the power method applied twice.
% Output: Smallest eigen value
% Date: January 29, 2018
% Supports only square matrices

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

out_file = fopen('/Users/kanikas/Documents/MatrixFree/Matrix-free/ConditionNo/eig_vals_cond_num_iter4_rel_err.csv','a') ;
fprintf(out_file,'Filename, Dimension , True lagest Eigen Value, Eigen Value from Power method, Absolute Error, Relative Error for largest Eigen Value, True smallest Eigen Value, Smallest eigen value from double power method, Relative error for smallest Eigen value, Condition No. from power method, True Condition number, Relative error for condition no., Symmetric? \n');


for i = 4: length(all_names)
    if all_files(i).bytes > 1000000 % less than 1 MB for now
        continue
    else 
        in_file = strcat(file_location,'/', all_files(i).name);
        A = load(in_file);
        sprintf('Computing values for file: %s', all_files(i).name)
        %Step 1: Get A
        A = A.Problem.A;
        if size(A,1)~= size(A,2) % only square matrices are handled
            continue
        else 
            n = size(A,1);
            w = rand(n,1); % choosing a random vector for initial guess
            w = w/norm(w);
            cntr = 0;
            symmetricity = issymmetric(A);

            %Step 2: Get largest eigen value of A by power method
            [lamda_max_A, true_eig_val_max_A] = power_method(A,w); 
            sprintf('True Largest Eigen value: %f', true_eig_val_max_A) % largest eigen value
            sprintf('Eigen value from power method: %f', lamda_max_A)
            abs_err = abs(true_eig_val_max_A) - abs(lamda_max_A);
            abs_err = abs(abs_err)
            rel_err = abs(abs_err / true_eig_val_max_A)

            %Step 3: Get B = A - (lamda_max * I) %shift
            I = speye(n); % for very large Identity matrix eye(n) doesnt work: Has error: Error using eye Requested 60000x60000 (26.8GB) array exceeds maximum array size preference. Creation of arrays greater than this limit may take a long time and cause MATLAB to become unresponsive. See array
            %size limit or preference panel for more information.
            B = A - (lamda_max_A * I); 

            %Step 4: Get lamda_max of B by Power method
            [lamda_max_B, true_eig_val_max_B] = power_method(B,w);

            %Step 5: Get smallest Eigen value lamda_min of A = (lamda_max of B) + (lamda_max of A) 
            lamda_min_A = lamda_max_B + lamda_max_A;

            try
                true_smallest_eigen_val = eigs(A, 1, 'sr');
            catch % for the case when eigs() fails to converge due to high condition number, reduce the tolerance of EIGS
                opts.tol = 1e-2; 
                true_smallest_eigen_val = eigs(A, 1, 'sr', opts);
            end 

            sprintf('True smallest eigen value of A: %f ', true_smallest_eigen_val)
            sprintf('Smallest Eigen value from double power method approach: %f ', lamda_min_A)
            abs_err_small = abs(true_smallest_eigen_val) - abs(lamda_min_A);
            abs_err_small = abs(abs_err_small)
            rel_err_small = abs(abs_err_small / true_smallest_eigen_val)
            condn_num = abs(lamda_max_A)/abs(lamda_min_A);
            sprintf('Condition number from double power method approach: %f', condn_num)
            true_condn_num = cond(A);
            abs_err_condn = abs(true_condn_num) - abs(condn_num);
            abs_err_condn = abs(abs_err_condn)
            sprintf('True Condition number: %f', true_condn_num)
            rel_err_cond = abs(abs_err_condn/true_condn_num);
            fprintf(out_file, '%s, %i, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f \n', all_files(i).name, n, true_eig_val_max_A, lamda_max_A, abs_err, rel_err, true_smallest_eigen_val, lamda_min_A, rel_err_small, condn_num, true_condn_num, rel_err_cond, symmetricity);
       end
    end
end

fclose(out_file);

function [lamda, true_eig_val] = power_method(A,w) 
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

    true_eig_val = val;
   
end


