load M1.txt;
M = spconvert(M1);
load b.txt;
load x1.txt;
x = pcg (M, b, 1.e-6);
error_x = abs(x - x1);
percentageDifference = error_x ./ x1;