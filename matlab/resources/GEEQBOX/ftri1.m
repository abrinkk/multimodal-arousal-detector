function fa = ftri1(x,z,n,ni,nstart);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Written by: Sarah Ratcliffe
% Based on Stata programs by Justine Shults
%
% COMMAND:  fa = ftri1(x,z,n,ni,nstart);
% ACTION:   Computes the stage 1 QLS estimates for the Tri-diagonal 
%           correlation structure assumption.
%
% SEE ALSO:  qls.m
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fa = 0;
for i=1:n
  nii = ni(i);
  zi = z(nstart(i):nstart(i)+nii-1,1);
  % Set Ri = tr-diagonal identity matrix;
  temp = ones(nii,1);
  Ri = spdiags([temp,x*temp,x*temp],[0 1 -1],nii,nii);  % Correlation matrix
  % Set dRi = tr-diagonal derivative matrix;
  dRi = spdiags([zeros(nii,1),temp,temp],[0 1 -1],nii,nii);  % Derivative matrix
  % Calculate z'*d(Rinv,alpha)*z;
  fa = fa - zi'*inv(Ri)*dRi*inv(Ri)*zi;
end;
