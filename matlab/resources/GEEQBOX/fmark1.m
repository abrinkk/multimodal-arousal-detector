function fa = fmark1(x,z,zlag,zzl,eij);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Written by: Sarah Ratcliffe
% Based on Stata programs by Justine Shults
%
% COMMAND:  fa = fmark1(x,z,zlag,zzl,eij,);
% ACTION:   Computes the stage 1 QLS estimates for the Markov correlation
%           structure assumption.
%
% SEE ALSO:  qls.m
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

aeij = x.^(eij);
a2eij = x.^(2*eij);
aeij1 = x.^(eij-1);
a2eij1 = x.^(2*eij-1);

fa1 = 2.*eij.*a2eij1.*(z.^2 - 2.*aeij.*zzl + zlag.^2) ./ ((1 - a2eij).^2);
fa2 = 2.*eij.*aeij1.*zzl ./ ((1-a2eij));
fa = nansum(fa1 - fa2);
