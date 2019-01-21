function [bhat,alpha,results] = qls(id,y,t,X,Family,CorrStruct,varnames,tol,maxit);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Written by: Sarah Ratcliffe
% Based on Stata programs by Justine Shults
% Copywrite 2006
%
% COMMAND:  [beta,alpha,results] = qls(id,y,t,X,Family,CorrStruct,varnames,tol,maxit);
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
%               1 | ar1:     AR(1) structure   
%               2 | markov:  Markov structure  (default)
%               3 | equi:    Equicorrelated structure
%               4 | tri:     Tri-diagonal structure
%		        5 | un:	     Unstructured correlation matrix
%		        6 | ind:	 Working independent correlation matrix
%    varnames = Optional cell array containing variable names associated
%               with the variables entered in X. (default = {'1','2',...})
%         tol = Convergence tolerance. (default = 0.00001)
%       maxit = Maximum number of iterations. (default = 100)
%
% OUTPUTS: bhat    = estimated effects parameters corresponding to X.
%          alpha   = estimated true correlation parameter.
%          results = structured array containing the estimates.
%               results.robust = estimates using robust covariance matrix.
%               results.model = estimates using model covariance matrix.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% NB: This program uses the following special m-files:
% fmark1.m  fmark2.m  fequi1.m
% counti.m

if nargin<9; maxit = 100; end;    % Maximum number of iterations;
if nargin<8; tol = 0.0000001; end;  % Convergence tolerance;
if nargin<6; CorrStruct=2; end;
if nargin<5; Family=1; end;
if nargin<4;
    disp('ERROR: At least 4 arguments must be specified: qls(id,y,t,X)');
    return;
end;
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

while (wtol>tol)&(iter<=maxit);
    
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

  % Generate correlation parameter estimate based on specified Correlation Structure;
  % Stage 1 estimates of alpha and correlation matrix, R;
  Z1 = [];
  Z2 = [];
  switch lower(CorrStruct);
  case {1,'ar1'};
    if iter==1; 
        disp('AR(1) Correlation structure assumed'); 
        disp(['Initial estimate of beta = [' num2str(bhat') ']']);
    end;
    a = nansum(z.^2 + zlag.^2);
    b = nansum(z.*zlag);
    c = nansum((z - zlag).^2);
    d = nansum((z + zlag).^2);
    alpha = (a - sqrt(c*d))/(2*b);  % Stage 1 estimate  
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
    end;
    
  case {2,'markov'};
    if iter==1; disp('Markov Correlation structure assumed'); 
        disp(['Initial estimate of beta = [' num2str(bhat') ']']);
    end;
    a = .00001;  % Starting value for lower interval limit.
    b = .99999;  % Starting value for upper interval limit.
    zzl = z.*zlag;   % zzl_ij = z_ij * z_ij-1
    alpha = fzero(@fmark1,[a b],[],z,zlag,zzl,eij);  % Stage 1 estimate
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
    end;

  case {3,'equi'};
    if iter==1; disp('Equicorrelated structure assumed'); 
                disp(['Initial estimate of beta = [' num2str(bhat') ']']);
    end;
    maxn = max(ni);
    a = -1/(maxn - 1) + 0.001;  % Lower limit for estimate
    b = .99999;                 % Upper limit for estimate
    p1 = sum(z.^2);
    for i=1:n-1;
        p2(i,1) = sum(z(nstart(i):nstart(i+1)-1));
    end;
    p2(n,1) = sum(z(nstart(n):totobs));
    alpha = fzero(@fequi1,[a b],[],p1,p2,ni);   % Stage 1 estimate
    for i=1:n;
       nii = ni(i);
       hui = huij(nstart(i):nstart(i)+nii-1,1);
       zi = z(nstart(i):nstart(i)+nii-1,1);
       Xi = X(nstart(i):nstart(i)+nii-1,:);
       Ri = alpha*ones(nii)+(1-alpha)*eye(nii);   % Correlation matrix
       SSi = chol(inv(Ri));
       Ai = diag(hui);
       Z1 = [Z1; SSi*zi];
       Z2 = [Z2; SSi*Ai*Xi];
    end;
    
  case {4,'tri'};
    if iter==1; 
        disp('Tri-diagonal Correlation structure assumed'); 
        disp(['Initial estimate of beta = [' num2str(bhat') ']']);
    end;
    maxn = max(ni);
    a = -1/(2*sin((maxn-1)/(maxn+1)*pi/2)) + 0.000001;  % Define lower value of feasible region for alpha.
    b = 1/(2*sin((maxn-1)/(maxn+1)*pi/2)) - 0.000001;  % Define upper value of feasible region for alpha.
    alpha = fzero(@ftri1,[a b],[],z,n,ni,nstart);   % Stage 1 estimate
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
    end;

  case {5,'un'};
    if iter==1; 
        disp('Unstructured Correlation matrix assumed'); 
        disp(['Initial estimate of beta = [' num2str(bhat') ']']);
    end;
    maxn = max(ni);
    phihat = sum(z.^2)/(totobs-p);
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
    end;
    
  case {6,'ind'};
    if iter==1; disp('Working Independent Correlation matrix assumed'); end;
    CorrStruct=6;
    alpha=0;
    
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
    wtol = 0;
  end;
%  format long g;
%  disp([iter; alpha; bhat])
%  format short;
  iter = iter+1;
end;

disp(' ');
disp(['Stage 1 estimate of alpha = ' num2str(alpha)]);
disp(['Stage 1 estimate of beta = [' num2str(bhat') ']']);

%%%%%%%%%%%%%%%%%%%%%%%%
% Stage 2 Estimation;  %
%%%%%%%%%%%%%%%%%%%%%%%%

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

Z1 = [];
Z2 = [];
SS = zeros(totobs,max(ni));

% Generate correlation parameter estimate based on specified Correlation Structure;
% Stage 2 esimate of phat
  switch lower(CorrStruct);
  case {1,'ar1'};
    phat = (2*alpha)/(1 + alpha^2);
    for i=1:n;
       nii = ni(i);
       ti = t(nstart(i):nstart(i)+nii-1,1);
       hui = huij(nstart(i):nstart(i)+nii-1,1);
       zi = z(nstart(i):nstart(i)+nii-1,1);
       Xi = X(nstart(i):nstart(i)+nii-1,:);
       Ri = phat.^(abs(kron([1:nii]',ones(1,nii)) - kron(ones(nii,1),[1:nii])));  % Correlation matrix
       SSi = chol(inv(Ri));
       Ai = diag(hui);
       Z1 = [Z1; SSi*zi];
       Z2 = [Z2; SSi*Ai*Xi];
       SS(nstart(i):nstart(i)+nii-1,1:nii) = SSi';
    end;

  case {2,'markov'};
    a = .00001;  % Starting value for lower interval limit.
    b = .99999;  % Starting value for upper interval limit.
    out1 = fmark2(a,alpha,eij);
    out2 = fmark2(b,alpha,eij);
    if sign(out1)~=sign(out2)
      phat = fzero(@fmark2,[a b],[],alpha,eij);
    else  
      disp('Warning: No stage 2 estimate possible as function to be solved has no zero over valid region');
      disp('Will use stage 1 estimate for alpha');
      phat = alpha;
    end  
    for i=1:n;
       nii = ni(i);
       ti = t(nstart(i):nstart(i)+nii-1,1);
       hui = huij(nstart(i):nstart(i)+nii-1,1);
       zi = z(nstart(i):nstart(i)+nii-1,1);
       Xi = X(nstart(i):nstart(i)+nii-1,:);
       Ri = phat.^(abs(kron(ti,ones(1,nii)) - kron(ones(nii,1),ti')));  % Correlation matrix
       SSi = chol(inv(Ri));
       Ai = diag(hui);
       Z1 = [Z1; SSi*zi];
       Z2 = [Z2; SSi*Ai*Xi];
       SS(nstart(i):nstart(i)+nii-1,1:nii) = SSi';
    end;
    
  case {3,'equi'};
    if sum(ni==n)==n;  % Equal number of measurements for all subjects
        phat = (alpha^2*(n-2) + 2*alpha)/(alpha^2*(n-1)+1);
    else;
        f1 = 0; f2 = 0;
        for i=1:n
            nii = ni(i);
            if nii>1;
                f1 = f1 + (alpha^2*(nii-1)*(nii-2)*nii + 2*alpha*(nii-1)*nii)/((1 + alpha*(nii-1))^2);
                f2 = f2 + (nii*(nii-1)*(alpha^2*(nii-1)+1))/((1 + alpha*(nii-1))^2);
            end;
        end;
        phat = f1 / f2;
    end;
    for i=1:n;
       nii = ni(i);
       hui = huij(nstart(i):nstart(i)+nii-1,1);
       zi = z(nstart(i):nstart(i)+nii-1,1);
       Xi = X(nstart(i):nstart(i)+nii-1,:);
       Ri = phat*ones(nii)+(1-phat)*eye(nii);   % Correlation matrix
       SSi = chol(inv(Ri));
       Ai = diag(hui);
       Z1 = [Z1; SSi*zi];
       Z2 = [Z2; SSi*Ai*Xi];
       SS(nstart(i):nstart(i)+nii-1,1:nii) = SSi';
    end;
    
  case  {4,'tri'};
    maxn = max(ni);
    a = -1/(2*sin((maxn-1)/(maxn+1)*pi/2)) + 0.000001;  % Define lower value of feasible region for alpha.
    b = 1/(2*sin((maxn-1)/(maxn+1)*pi/2)) - 0.000001;  % Define upper value of feasible region for alpha.
    out1 = sign(ftri2(a,alpha,z,n,ni,nstart));
    out2 = sign(ftri2(b,alpha,z,n,ni,nstart));
    if out1~=out2
      phat = fzero(@ftri2,[a b],[],alpha,z,n,ni,nstart);   % Stage 2 estimate
    else  
      phat = fminsearch(@ftri2a,(a+b)/2,[],alpha,z,n,ni,nstart);   % Stage 2 estimate
      disp(' ');
      disp('Warning: Exact stage 2 estimate unknown as function to be solved has no zero over valid region');
      disp('Alpha value that results in function value closest to zero within valid region will be used.');
      if phat>b; phat=b; end;
      if phat<a; phat=a; end;
    end  
    for i=1:n;
       nii = ni(i);
       hui = huij(nstart(i):nstart(i)+nii-1,1);
       zi = z(nstart(i):nstart(i)+nii-1,1);
       Xi = X(nstart(i):nstart(i)+nii-1,:);
       temp = ones(nii,1);
       Ri = spdiags([temp,phat*temp,phat*temp],[0 1 -1],nii,nii);  % Correlation matrix
       SSi = chol(inv(Ri));
       Ai = diag(hui);
       Z1 = [Z1; SSi*zi];
       Z2 = [Z2; SSi*Ai*Xi];
       SS(nstart(i):nstart(i)+nii-1,1:nii) = SSi';
    end;

  case {5,'un'};
    maxn = max(ni);
    phihat = sum(z.^2)/(totobs-p);
    phat = 0;
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
    alpha=0;
    
  otherwise
    disp('ERROR: Invalid correlation structure specified');
    return;
  end;

if CorrStruct~=6;  
  gammahat = regress(Z1,Z2);
  bhat = bhat + gammahat;
  alpha = phat;
end;

disp(' ');
disp(['Stage 2 estimate of alpha = ' num2str(alpha)]);
disp(['Stage 2 estimate of beta = [' num2str(bhat') ']']);
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

disp(['QLS estimate of scale parameter = ' num2str(tau)]);
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
