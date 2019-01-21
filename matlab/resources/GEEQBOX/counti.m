function [ni,n] = counti(x);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% copywrite Sarah Ratcliffe 2000
%
% COMMAND:  [ni,n] = counti(x)
% ACTION:  For a vector of non-unique elements, x,
%          returns the number of times each unique
%          element is contained in x in ni. The
%          number of unique elements is optionally
%          returned in n.
%    Note: the non-unique elements of x must be 
%          grouped together and ni will contain in
%          numbers of each non-unique element in
%          the order in which they appear in x.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if nargin<1
  disp('ERROR: missing input variable');
  return;
end;

[n,p] = size(x);
if (n>1) & (p>1)
  disp('ERROR: x must be a vector');
  return;
end;

%% Initial Values
p = n*p;	% Number of elements in x
ni(1) = 1;	% First subject has one element
n = 1;		% One unique subject to start
xi = x(1);	% First unique subject = first element in x

for i=2:p    % cycle through elements counting as we go
  if x(i)==xi;
    ni(n) = ni(n)+1;
  else
    n = n+1;
    ni(n) = 1;
    xi = x(i);
  end
end;


