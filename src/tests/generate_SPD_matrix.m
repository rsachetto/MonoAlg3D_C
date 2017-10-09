function generate_SPD_matrix(n, density)

A = sprandsym(n,density); % generate a random n x n matrix
A = A + n*speye(n);

X = rand(n, 1);

B = A*X;

save A.txt A
save B.txt B
save X.txt X

endfunction