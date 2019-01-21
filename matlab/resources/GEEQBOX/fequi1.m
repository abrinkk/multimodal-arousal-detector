 function fa = fequi1(x,p1,p2,p3);
 %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Written by: Sarah Ratcliffe
% Based on Stata programs by Justine Shults
%
% COMMAND:  fa = fequi1(x,p1,p2,p3);
% ACTION:   Computes the stage 1 QLS estimates for the equicorrelated
%           structure assumption. 
%
% INPUTS:    x = alpha to be solved for
%           p1 = sum_(ij) z_(ij)^2  (constant)
%           p2 = sum_j z_(ij)^2   (n x 1 vector)
%           p3 = n_i              (n x 1 vector)
%
% SEE ALSO:  qls.m
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

f1 = (p2.^2).*(1 + x^2*(p3-1));
f2 = (1 + x*(p3-1)).^2;
fa = p1 - sum(f1./f2);
