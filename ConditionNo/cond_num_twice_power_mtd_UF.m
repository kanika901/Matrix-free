% Computes the smallest eigen value with the eig() and the power method applied twice.
% Output: Smallest eigen value (smallest absolute value) for UFlorida matrices
% Date: January 29, 2018
% Supports only square matrices

% Pseudocode
% This approach uses Power Method twice
%1. Get A
%2. Get lamda_max of A by Power method
%3. B = A - (lamda_max * I) %shift
%4. Get lamda_max of B by Power method
%5. lamda_min of A = (lamda_max of B) + (lamda_max of A) 

cd '/Users/kanikas/Documents/mat_files_1014'
%file_location = '/Users/kanikas/Documents/research/data_(mat)_matlab_format/mat_files';
file_location = '/Users/kanikas/Documents/mat_files_1014';
all_files = dir(file_location);
all_names = { all_files.name };
global in_file;

out_file = fopen('/Users/kanikas/Documents/MatrixFree/Matrix-free/ConditionNo/eig_vals_cond_num_iter4_rel_err_all_matr_final_v2.csv','a') ;
fprintf(out_file,'Filename, Dimension, True largest Eigen Value, Eigen Value from Power method, Absolute Error, Relative Error for largest Eigen Value, True smallest Eigen Value, Smallest eigen value from double power method, Relative error for smallest Eigen value, Condn No. from lamda_max/lamda_min, True Condn no. (condest(A 2)), condest(A 1), cond(A 2), cond(A 1), Relative error for condition no., Symmetric? \n');

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

            %true_smallest_eigen_val = eigs(A, 1, 'sr');
            %We want smallest magnitude eigen value
            true_smallest_eigen_val = min(abs(eigs(A))); %eigs(A, 1, 'SM'); 
            
            sprintf('True smallest eigen value of A: %f ', true_smallest_eigen_val)
            sprintf('Smallest Eigen value from double power method approach: %f ', lamda_min_A)
            abs_err_small = abs(true_smallest_eigen_val) - abs(lamda_min_A);
            abs_err_small = abs(abs_err_small)
            rel_err_small = abs(abs_err_small / true_smallest_eigen_val)
            if lamda_min_A == 0
                condn_num = Inf;
            else
                condn_num = abs(lamda_max_A)/abs(lamda_min_A);
            end
            sprintf('Condition number from lamda_max/lamda_min: %f', condn_num)
            true_condn_num = condest(A,2);
            abs_err_condn = abs(true_condn_num) - abs(condn_num);
            abs_err_condn = abs(abs_err_condn)
            sprintf('True Condition number: %f', true_condn_num)
            rel_err_cond = abs(abs_err_condn/true_condn_num);
            cond_est_1 = condest(A,1);
            cond_2 = cond(A,2);
            cond_1 = cond(A,1);
            fprintf(out_file, '%s, %i, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f \n', all_files(i).name, n, true_eig_val_max_A, lamda_max_A, abs_err, rel_err, true_smallest_eigen_val, lamda_min_A, rel_err_small, condn_num, true_condn_num, cond_est_1, cond_2, cond_1, rel_err_cond, symmetricity);
       end
    end
movefile(in_file,'../DONE')
end
fclose(out_file);

function [lamda, true_eig_val] = power_method(A,w) 
    global in_file
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
    
    try % to handle the error --> with ARPACK routine dneupd: dnaupd did not find any eigenvalues to sufficient accuracy.
         val = eigs(A, 1, 'LM'); %max(abs(eigs(A)))
    catch
         opts.tol = 1e-6;
         try
             val = eigs(A, 1, 'LM', opts);
         catch
             movefile(in_file,'../notUsed')
         end 
    end
%     val_with_sign = max(eigs(A));
%     if val == val_with_sign
%        val_sign = 1; % positive
%     else   
%        val_sign = -1; % negative      
%     end
%        val = val * val_sign;
    true_eig_val = val;  
end


