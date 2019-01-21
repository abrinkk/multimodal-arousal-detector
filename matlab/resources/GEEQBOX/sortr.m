function [x,Z] = sortr(x,Zold);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% copyright Sarah Ratcliffe 2001
%
% COMMAND:  [x,Z] = sortr(x,Z);
% ACTION:  Sorts repeated measures data (Z) according
%          to the order needed to sort the subject 
%          unique elements of the vector x.
%          The ordering of x must correspond to the
%          ordering of the subject repeated measurements
%          in Z, with the first column of Z being the
%          unique subject ID. 
%          Main use with rpat.m, sorting data according
%          to pattern ordering.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[ni,n] = counti(Zold(:,1));

if length(x) ~= n
  disp('ERROR: more subjects in Z than elements in x');
  return;
end;

[x,b] = sort(x);
Z = [];

for i=1:n
  ref = b(i);
  nis = sum(ni(1:ref-1))+1;
  nie = sum(ni(1:ref));
  Z = [Z; Zold(nis:nie,:)];
end;

  
  
