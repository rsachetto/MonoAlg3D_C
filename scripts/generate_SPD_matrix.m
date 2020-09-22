function generate_SPD_matrix(n, density)

A = sprandsym(n,density); % generate a random n x n matrix
A = A + n*speye(n);

X = rand(n, 1);

B = A*X;

save A4.txt A
save B4.txt B
save X4.txt X

endfunction
