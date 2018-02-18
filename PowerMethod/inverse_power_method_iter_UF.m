% Computes the smallest eigen value with the eig() and the inverse power method.
% Output: Smallest eigen value
% Date: January 22, 2018

file_location = '/Users/kanikas/Documents/research/data_(mat)_matlab_format/mat_files';
all_files = dir(file_location);
all_names = { all_files.name };

out_file = fopen('/Users/kanikas/Documents/MatrixFree/Matrix-free/PowerMethod/eig_val_smallest_results_iter4_rel_err.csv','a') ;
fprintf(out_file,'Filename, Dimension , No. of Iterations, Eigen Value, Eigen Value from Power method, Absolute Error, Relative Error \n');

for i = 4: length(all_names)
    if all_files(i).bytes > 1000000 % less than 1 MB for now
        continue
    else 
        in_file = strcat(file_location,'/', all_files(i).name);
        A = load(in_file);
        sprintf('Computing values for file: %s', all_files(i).name)
        A = A.Problem.A;
        n = size(A,1);
        y = rand(n,1); % choosing a random vector for initial guess
        cntr = 0;
        iter = 1;
        %abs_err = 0;
        %rel_err = 0;
        if size(A,1)~= size(A,2) % only square matrices are handled
            continue
        else
            b = y;
            while (iter < 5) %fixed # iterations
                cntr = cntr + 1;
                iter = iter + 1;
                w = A\b;
                norm_w = norm(w);
                lamda = 1/norm_w;
                b = w/norm_w;
            end

            sprintf('Dimension: %i', n)
            try
                true_eig_val = eigs(A,1,'sm');
                sprintf('Eigen value: %f', true_eig_val) % smallest eigen value
                sprintf('Eigen value from inverse power method: %f', lamda)
                abs_err = abs(true_eig_val) - abs(lamda);
                abs_err = abs(abs_err);
                rel_err = abs(abs_err / true_eig_val);
                fprintf(out_file, '%s, %i, %i, %f, %f, %f, %f \n', all_files(i).name, n, iter - 1, true_eig_val, lamda, abs_err, rel_err);
            end
            %fprintf(out_file, '%s, %i, %i, %i, %i, %i \n', all_files(i).name, n, iter, true_eig_val, lamda, abs(true_eig_val - lamda));
    end
    end
end
fclose(out_file);