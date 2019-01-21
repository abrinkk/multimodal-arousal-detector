function fa = ftri2(x,alpha,z,n,ni,nstart);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Written by: Sarah Ratcliffe
% Based on Stata programs by Justine Shults
%
% COMMAND:  fa = ftri2(x,alpha,z,n,ni,nstart,signind);
% ACTION:   Computes the stage 2 QLS estimates for the Tri-diagonal 
%           correlation structure assumption.
%
% SEE ALSO:  qls.m
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fa = 0;
for i=1:n
  nii = ni(i);
  zi = z(nstart(i):nstart(i)+nii-1,1);
  % Set Ri = tr-diagonal matrix at stage 1 estimate;
  temp = ones(nii,1);
  Ri = spdiags([temp,alpha*temp,alpha*temp],[0 1 -1],nii,nii);  % Correlation matrix
  % Set dRi = tr-diagonal derivative matrix;
  dRi = spdiags([zeros(nii,1),temp,temp],[0 1 -1],nii,nii);  % Derivative matrix
  % Set Rinew = tr-diagonal matrix at stage 2 estimate;
  temp = ones(nii,1);
  Rinew = spdiags([temp,x*temp,x*temp],[0 1 -1],nii,nii);  % Correlation matrix
  % Calculate trace(d(Rinv,alpha)*R(new));
  fa = fa + trace(inv(Ri)*dRi*inv(Ri)*Rinew);
end;
