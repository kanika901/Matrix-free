% Computes the larget eigen value with the eig() and the power method.
% Output: Largest eigen value with power method, testing UF matrices for accuracy evaluation
% Date: December 1, 2017

file_location = '/Users/kanikas/Documents/research/data_(mat)_matlab_format/mat_files';
all_files = dir(file_location);
all_names = { all_files.name };

out_file = fopen('/Users/kanikas/Documents/MatrixFree/Matrix-free/PowerMethod/eig_val_results.csv','a') ;
fprintf(out_file,'Filename, Dimension , No. of Iterations, Eigen Value, Eigen Value from Power method, Error \n');

for i = 4: length(all_names)
    if all_files(i).bytes > 1000000 % less than 1 MB for now
        continue
    else 
        in_file = strcat(file_location,'/', all_files(i).name);
        A = load(in_file);
        sprintf('Computing values for file: %s', all_files(i).name)
        A = A.Problem.A;
        n = size(A,1);
        w = rand(n,1); % choosing a random vector for initial guess
        w = w/norm(w);
        cntr = 0;
        iter = 0;
        if size(A,1)~= size(A,2) % only square matrices are handled
            continue
        else
            while (iter < 3) %fixed # iterations
                cntr = cntr + 1;
                iter = iter + 1;
                b = A * w; %  matrix-vector product 
                lamda = max(abs(b));
                w = b/lamda; %x1 = x0/?
            end
            val = eigs(A);
            sprintf('Dimension: %i', n)
            sprintf('Eigen value: %f', val(1)) % largest eigen value
            sprintf('Eigen value from power method: %f', lamda)
            sprintf('No. of matrix-vector products: %i', cntr) 
            fprintf(out_file, '%s, %i, %i, %i, %i, %i \n', all_files(i).name, n, iter, val(1), lamda, abs(val(1) - lamda));
    end
    end
end
fclose(out_file);