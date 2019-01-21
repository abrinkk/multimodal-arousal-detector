function plotsub(id,x,y,S);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Line plot of x,y data for individual subjects
% Written by: Sarah Ratcliffe
%
% COMMAND:  plotsub(id,x,y,S);
%
% ACTION: same as the plot command except for separate
%         line plots will be produced for each subject.
%         
% INPUTS: id = (nx1) vector containing the unique subject ids. This
%              must have all data for a subject grouped together.
%          x = (nx1) vector containing x variable for plot.
%          y = (nx1) vector containing y variable for plot.
%          S = optional string argument to be passed to plot to 
%              control plotting options.
%
% SEE ALSO: plot.m
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[ni,n] = counti(id);

nend = cumsum(ni);
nstart = [1 nend+1];

if nargin==3;
  S='';
end;
plot(x(nstart(1):nend(1)),y(nstart(1):nend(1)),S);
hold on;
for i=2:n
  plot(x(nstart(i):nend(i)),y(nstart(i):nend(i)),S);
end;