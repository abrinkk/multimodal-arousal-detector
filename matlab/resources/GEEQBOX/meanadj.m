function xm = meanadj(x,groups);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Written by: Sarah Ratcliffe
% Last update: 12 Sep 1997
%
% 	function xm = meanadj(x,groups);
%
% Action: mean adjusts the data matrix x for each of the groups individually.
%
% Input:       x = matrix of data values. (longitudinal data)
%		    A subject's values go across the rows.
%	  groups = a vector containing the number of subjects in each group.
%
% Output:     xm = a matrix of mean adjusted data values.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if nargin == 2;
ng = length(groups);
startg = 1;
endg = 0;

for i=1:ng;
  ni = groups(i);
  endg = endg + ni;

  xg = x(startg:endg,:);		% extract group data
  meang = mean(xg);			% calculate mean of group

  xg = xg-(ones(ni,1)*meang);		% adjust each subject
  
  xm(startg:endg,:) = xg;
  startg = startg + ni;
end

else
disp('Must have two arguements - x and groups');
end;
