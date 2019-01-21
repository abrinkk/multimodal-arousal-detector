function fa = fmark2(x,alpha,eij);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Written by: Sarah Ratcliffe
% Based on Stata programs by Justine Shults
%
% COMMAND:  fa = fmark2(x,alpha,eij,);
% ACTION:   Computes the stage 2 QLS estimates for the Markov correlation
%           structure assumption.
%
% SEE ALSO:  qls.m
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

alpeij = alpha.^(eij);
aeij = x.^(eij);

fa1 = 2.*eij.*(alpeij.^2) - aeij.*eij.*(alpeij + (alpeij.^3));
fa2 = (1 - (alpeij.^2)).^2;
fa = nansum(fa1 ./ fa2);
