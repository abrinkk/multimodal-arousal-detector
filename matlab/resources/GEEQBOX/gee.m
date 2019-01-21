function [bhat,alpha,results] = gee(id,y,t,X,Family,CorrStruct,varnames,tol,maxit);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Written by: Sarah Ratcliffe
% Based on Stata programs by Justine Shults
% Copywrite 2006
%
% COMMAND:  [beta,alpha,results] = gee(id,y,t,X,Family,CorrStruct,varnames,tol,maxit);
% ACTION:   Computes Quassi Least squares (QLS) estimates for repeated
%           measures data y, measured at times t, on subjects id, regressed
%           on covariates X, assuming data is from a Family distribution 
%           with CorrStruct type of correlation.
%
% INPUTS:  id = Column of unique subject identifiers used to indicate
%               measurements from the same subjects.
%           y = Column of outcome variable measurements corresponding to id.
%           t = Column of measurement times corresponding to y.
%               NB: AT THE MOMENT ASSUMING EQUISPACED DATA.
%           X = Matrix of covariates corresponding to y and t. This must 
%               contain a column of 1's in order to include a constant term 
%               in the model.
%          Family = Assumed distribution of the data.
%               1 | n:       Normal distribution  (default)
%               2 | b:       Bernoulli distribution
%               3 | p:       Poisson
%          CorrStruct = Assumed correlation structure
%               1 | ar1:     AR(1) structure  (default) 
%               2 | markov:  Markov structure  
%               3 | equi:    Equicorrelated structure
%               4 | tri:     Tri-diagonal structure
%		        5 | un:	     Unstructured correlation matrix
%		        6 | ind:	 Working independent correlation matrix
%    varnames = Optional cell array containing variable names associated
%               with the variables entered in X. (default = {'1','2',...})
%         tol = Convergence tolerance. (default = 0.00001)
%       maxit = Maximum number of iterations. (default = 100)
%
% OUTPUTS:   bhat  = estimated effects parameters corresponding to X.
%            alpha = estimated true correlation parameters.
%          results = structured array containing the estimates.
%               results.robust = estimates using robust covariance matrix.
%               results.model = estimates using model covariance matrix.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% NB: This program uses the following special m-files:
% fmark1.m  fmark2.m  fequi1.m
% counti.m

if nargin<9; maxit = 100; end;    % Maximum number of iterations;
if nargin<8; tol = 0.00001; end;  % Convergence tolerance;
if nargin<6; CorrStruct=1; end;
if nargin<5; Family=1; end;

warning off MATLAB:divideByZero
warning off MATLAB:fzero:UndeterminedSyntax

% Check dimensions of input variables.
[totobs,ct] = size(id);
if totobs==1;  % Change data from row to column vector;
  id = id';
  totobs=ct;
  ct=1;
end;
if ct~=1  % Not 1 column
  disp('ERROR: id must contain a single column of data listing subject ids');
  return;
end;
[rt,ct] = size(y);
if rt==1;  % Change data from row to column vector;
  y = y';
  rt=ct;
  ct=1;
end;
if ct~=1 | rt~=totobs % Not 1 column of same length as id
  disp('ERROR: y must contain a single column of data from the outcome variable');
  disp('       and must be the same length as id');
  return;
end;
[rt,ct] = size(t);
if rt==1;  % Change data from row to column vector;
  t = t';
  rt=ct;
  ct=1;
end;
if ct~=1 | rt~=totobs  % Not 1 column of same length as id
  disp('ERROR: t must contain a single column of time data of same lenght as y');
  return;
end;
[rt,ct] = size(X);
if rt~=totobs;
  disp('ERROR: X must have same number of rows as length of y');
  return;
end;
[totobs,p] = size(X);

% Determine number of observations per subject.
[ni,n] = counti(id);
ni = ni';
nstart = [1; cumsum(ni)+1];
nstart = nstart(1:n);

% Generate starting values for beta = bhat;
switch lower(Family);
case {1,'n'}
    disp('Normal distribution family assumed');
    bhat = glmfit(X,y,'normal','identity','on',[],[],'off');	
case {2,'b'};
    disp('Bernoulli distribution family assumed');
    bhat = glmfit(X,[y ones(totobs,1)],'binomial','logit','on',[],[],'off');
case {3,'p'};
    disp('Poisson distribution family assumed');
    bhat = glmfit(X,y,'poisson','log','on',[],[],'off');
otherwise
  disp('ERROR: Invalid distribution family specified');
  return;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start Stage 1 estimation loop; %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wtol = 1;
iter = 1;

tlag = [NaN; t(1:totobs-1)];
tlag(nstart) = NaN;
eij = t - tlag;  % e_ij = t_ij - t_ij-1

SS = zeros(max(ni),max(ni),n);
K1p = (totobs-n)-p;

while (wtol>tol)&(iter<=maxit);
%disp(['GEE estimate of beta = [' num2str(bhat') ']']);

  xtb = X*bhat;
  % Create z_ij (Pearson residuals), h(u_ij) vectors
  switch lower(Family);
  case {1,'n'}
    z = y - xtb;
    huij = ones(totobs,1);
  case {2,'b'};
    uij = exp(xtb)./(1+exp(xtb));
    huij = sqrt(uij.*(1-uij));
    z = (y - uij)./huij;
  case {3,'p'};
    uij = exp(xtb);
    huij = sqrt(uij);
    z = (y - uij)./huij;
  otherwise
  end;

  zlag = [NaN; z(1:totobs-1,1)];
  zlag(nstart) = NaN;

  % Estimate dispersion paramater, phi;
  phihat = sum(z.^2)/(totobs-p);
  
  % Generate correlation parameter estimate based on specified Correlation Structure;
  % Stage 1 estimates of alpha and correlation matrix, R;
  Z1 = [];
  Z2 = [];
  SS = zeros(totobs,max(ni));

  switch lower(CorrStruct);
  case {1,'ar1'};
    if iter==1; disp('AR(1) Correlation structure assumed'); end;
    
    alpha = nansum(z.*zlag)/(K1p*phihat);

    for i=1:n;
       nii = ni(i);
       hui = huij(nstart(i):nstart(i)+nii-1,1);
       zi = z(nstart(i):nstart(i)+nii-1,1);
       Xi = X(nstart(i):nstart(i)+nii-1,:);
       Ri = alpha.^(abs(kron([1:nii]',ones(1,nii)) - kron(ones(nii,1),[1:nii])));  % Correlation matrix
       SSi = chol(inv(Ri));
       Ai = diag(hui);
       Z1 = [Z1; SSi*zi];
       Z2 = [Z2; SSi*Ai*Xi];
       SS(nstart(i):nstart(i)+nii-1,1:nii) = SSi';
    end;
    
  case {2,'markov'};
    if iter==1; disp('Markov Correlation structure assumed'); end;
    
    alpha = nansum(z.*zlag)/(K1p*phihat);

    for i=1:n;
       nii = ni(i);
       ti = t(nstart(i):nstart(i)+nii-1,1);
       hui = huij(nstart(i):nstart(i)+nii-1,1);
       zi = z(nstart(i):nstart(i)+nii-1,1);
       Xi = X(nstart(i):nstart(i)+nii-1,:);
       Ri = alpha.^(abs(kron(ti,ones(1,nii)) - kron(ones(nii,1),ti')));  % Correlation matrix
       SSi = chol(inv(Ri));
       Ai = diag(hui);
       Z1 = [Z1; SSi*zi];
       Z2 = [Z2; SSi*Ai*Xi];
       SS(nstart(i):nstart(i)+nii-1,1:nii) = SSi';
    end;

  case {3,'equi'};
    if iter==1; disp('Equicorrelated structure assumed'); end;
    maxn = max(ni);
    
    temp1 = 0;
    temp2 = 0;
    for i=1:n;
       nii = ni(i);
       zi = z(nstart(i):nstart(i)+nii-1,1);
       temp = zi * zi';
       temp1 = temp1 + sum(sum(temp)) - sum(diag(temp));
       temp2 = temp2 + nii*(nii-1);
    end;
    alpha = temp1 / (phihat*(temp2-p));

    for i=1:n;
       nii = ni(i);
       hui = huij(nstart(i):nstart(i)+nii-1,1);
       zi = z(nstart(i):nstart(i)+nii-1,1);
       Xi = X(nstart(i):nstart(i)+nii-1,:);
       Ri = alpha*ones(nii,nii)+(1-alpha)*eye(nii);   % Correlation matrix
       SSi = chol(inv(Ri));
       Ai = diag(hui);
       Z1 = [Z1; SSi*zi];
       Z2 = [Z2; SSi*Ai*Xi];
       SS(nstart(i):nstart(i)+nii-1,1:nii) = SSi';
    end;

  case {4,'tri'};
    if iter==1; disp('Tri-diagonal correlation structure assumed'); end;

    alpha = nansum(z.*zlag)/(K1p*phihat);
    for i=1:n;
       nii = ni(i);
       hui = huij(nstart(i):nstart(i)+nii-1,1);
       zi = z(nstart(i):nstart(i)+nii-1,1);
       Xi = X(nstart(i):nstart(i)+nii-1,:);
       temp = ones(nii,1);
       Ri = spdiags([temp,alpha*temp,alpha*temp],[0 1 -1],nii,nii);  % Correlation matrix
       SSi = chol(inv(Ri));
       Ai = diag(hui);
       Z1 = [Z1; SSi*zi];
       Z2 = [Z2; SSi*Ai*Xi];
       SS(nstart(i):nstart(i)+nii-1,1:nii) = SSi';
    end;
    
  case {5,'un'};
    if iter==1; disp('Unstructured Correlation matrix assumed'); end;
    maxn = max(ni);
    alpha = 0;
    rjk = zeros(maxn,maxn);
    for i=1:n;
       nii = ni(i);
       zi = z(nstart(i):nstart(i)+nii-1,1);
       rjk(1:nii,1:nii) = rjk(1:nii,1:nii) + zi*zi';
   end; 
   rjk = rjk/(phihat*(n-p));
   Rall = rjk+diag(1-diag(rjk));   % Correlation matrix
   for i=1:n;
       nii = ni(i);
       hui = huij(nstart(i):nstart(i)+nii-1,1);
       zi = z(nstart(i):nstart(i)+nii-1,1);
       Xi = X(nstart(i):nstart(i)+nii-1,:);
       Ri = Rall(1:nii,1:nii);
       SSi = chol(inv(Ri));
       Ai = diag(hui);
       Z1 = [Z1; SSi*zi];
       Z2 = [Z2; SSi*Ai*Xi];
       SS(nstart(i):nstart(i)+nii-1,1:nii) = SSi';
    end;
   
  case {6,'ind'};
    if iter==1; disp('Working Independent Correlation matrix assumed'); end;
    CorrStruct=6;
    alpha=1;
    
  otherwise
    disp('ERROR: Invalid correlation structure specified');
    return;
  end;

  % Estimate regression coefficients;
  if CorrStruct~=6;
    gammahat = regress(Z1,Z2);
    wtol = gammahat'*gammahat;
    bhat = bhat + gammahat;
  else
    wtol=0;    
  end;

  iter = iter+1;
end;

disp(' ');
disp(['GEE estimate of alpha = ' num2str(alpha)]);
disp(['GEE estimate of beta = [' num2str(bhat') ']']);
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%
% Other Stuff          %
%%%%%%%%%%%%%%%%%%%%%%%%

% Estimate of scalar parameter, phi which will be stored in the variable tau;
switch lower(Family);
case {1,'n'}
    sum1 = 0;
    sum2 = 0;
    for i=1:n;
       nii = ni(i);
       SSi = SS(nstart(i):nstart(i)+nii-1,1:nii);
       zi = z(nstart(i):nstart(i)+nii-1,1);
       temp = SSi*zi;
       sum1 = sum1 + 1/nii*temp'*temp;
       sum2 = sum2 + 1/nii*zi'*zi;
    end;
	tau = min(sum1,sum2)/n;
    
case {2,'b'};
    tau = 1;
case {3,'p'};
    tau = 1;
otherwise
end;

disp(['GEE estimate of scale parameter = ' num2str(tau)]);
disp(' ');

% P-values, 95% confidence intervals;
nwa = zeros(p,p);
newb = zeros(p,p);
for i=1:n
    nii = ni(i);
    Z1i = Z1(nstart(i):nstart(i)+nii-1,:);
    Z2i = Z2(nstart(i):nstart(i)+nii-1,:);
    nwa = nwa + Z2i'*Z2i;
    temp = Z2i'*Z1i;
    newb = newb + temp*temp';
end;
temp = inv(nwa);
robcov = temp*newb*temp;    % Robust covariance matrix;
modcov = tau*temp;          % Model based covariance matrix;

stderr = sqrt(diag(robcov));
cim = norminv(.975)*stderr;
ciup = bhat + cim;
cilo = bhat - cim;
zhat = stderr.^(-1).*bhat;
pval = 2*(1-normcdf(abs(zhat)));

cols = {'','','','','','95% CI','';'Variable','Beta','Std.Error','z value','p-value','low lim','up lim'};
if nargin<7; varnames= num2cell(1:p); end;

resr = [cols; [varnames' num2cell([bhat stderr zhat pval cilo ciup])]];
disp('Estimates based on ROBUST covariance matrix');
disp(resr)
results.robust = resr;

stderr = sqrt(diag(modcov));
cim = norminv(.975)*stderr;
ciup = bhat + cim;
cilo = bhat - cim;
zhat = stderr.^(-1).*bhat;
pval = 2*(1-normcdf(abs(zhat)));

resr = [cols; [varnames' num2cell([bhat stderr zhat pval cilo ciup])]];
disp('Estimates based on MODEL based covariance matrix');
disp(resr);
results.model = resr;
