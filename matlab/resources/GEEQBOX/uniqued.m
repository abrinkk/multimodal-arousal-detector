function list = uniqued(vec, d);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Written by: Sarah Ratcliffe
% Last Updated: 3 July 1997
%
% COMMAND  :  function list = uniqued(vec, d); 
%  ACTION  :  Returns the unique conponents of a vector to d
%		decimal places.
%
%   INPUT  :  vec = vector, that may contain repeats.
%		d = number of decimal places.
%
%  OUTPUT  :  list = vector of unique components.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

x = sort(vec);
x = x.*(10^d);
x = round(x);
x = x./(10^d);
list = unique(x);
