function kappa = fleisskappa(X)
%FLEISSKAPPA calculates Fleiss' Kappa measure.
%   kappa = FLEISSKAPPA(X) calculates the Fleiss' Kappa
%   statistic between vectors in X
%
%   Input:  X, inputs signals
%   Output: kappa, Fleiss' Kappa

% Number of samples
N = size(X,1);
% Number of raters
n = size(X,2);
% Different categories
C = unique(X(:));
% Number of categories
k = length(C);
% Number of raters who assigned the i-th subject to the j-th category
n_ij = arrayfun(@(x) sum(X == x,2), C, 'Un',0);
n_ij = [n_ij{:}];
% Proportion of all assignments which were to the j-th category
p_j = 1/(N*n)*sum(n_ij);
% Extent to which raters agree for the i-th subject
P_i = 1/(n*(n-1))*(sum(n_ij.^2,2) - n);
% The mean of P_i
P = 1/N*sum(P_i);
% Hypothetical probability of chance agreement
P_e = sum(p_j.^2);
% Fleiss' Kappa
kappa = (P - P_e)/(1 - P_e);
end