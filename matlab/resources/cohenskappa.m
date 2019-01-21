function [kappa,n_k1,n_k2,N] = cohenskappa(x1,x2)
%COHENSKAPPA calculates Cohen's Kappa measure.
%   [kappa,n_k1,n_k2,N] = COHENSKAPPA(x1,x2) calculates the Cohen's Kappa
%   statistic between vector x1 and x2.
%
%   Input:  x1, input signal 1
%           x2, input signal 2
%   Output: kappa, Cohen's Kappa
%           n_k1, number of samples in different category for signal 1
%           n_k2, number of samples in different category for signal 2
%           N, number of samples

% Different categories
C = unique([x1(:); x2(:)]);
% Number of samples
N = length(x1);
% Number of different categories for each scorer
n_k1 = arrayfun(@(x) sum(x1 == x), C);
n_k2 = arrayfun(@(x) sum(x2 == x), C);
% Hypothetical probability of chance agreement
p_e = 1/N^2*sum(n_k1.*n_k2);
% Accuracy
p_o = sum(x1 == x2)/N;
% Cohen's kappa
kappa = (p_o - p_e)/(1 - p_e);

end
