function [mse] = calc_mse(X, Xapp)
    D = abs(X-Xapp).^2;
    mse = sum(D(:))/numel(X);
