function [B,SE,PVAL,in,stats,nextstep,history] = stepwisefit_gee(allx,y,id,t,approach,varargin)
%STEPWISEFIT Fit regression model using stepwise regression
%   B=STEPWISEFIT(X,Y) uses stepwise regression to model the response variable
%   Y as a function of the predictor variables represented by the columns
%   of the matrix X.  The result B is a vector of estimated coefficient values
%   for all columns of X.  The B value for a column not included in the final
%   model is the coefficient that would be obtained by adding that column to
%   the model.  STEPWISEFIT automatically includes a constant term in all
%   models.
%
%   [B,SE,PVAL,INMODEL,STATS,NEXTSTEP,HISTORY]=STEPWISEFIT(...) returns additional
%   results.  SE is a vector of standard errors for B.  PVAL is a vector of
%   p-values for testing if B is 0.  INMODEL is a logical vector indicating
%   which predictors are in the final model.  STATS is a structure containing
%   additional statistics.  NEXTSTEP is the recommended next step -- either
%   the index of the next predictor to move in or out, or 0 if no further
%   steps are recommended.  HISTORY is a structure containing information
%   about the history of steps taken.
%
%   [...]=STEPWISEFIT(X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies one or
%   more of the following name/value pairs:
%
%     'inmodel'  A logical vector, or a list of column numbers, indicating which
%                predictors to include in the initial fit (default none)
%     'penter'   Max p-value for a predictor to be added (default 0.05)
%     'premove'  Min p-value for a predictor to be removed (default 0.10)
%     'display'  Either 'on' (default) to display information about each
%                step or 'off' to omit the display
%     'maxiter'  Maximum number of steps to take (default is no maximum)
%     'keep'     A logical vector, or a list of column numbers, indicating which
%                predictors to keep in their initial state (default none)
%     'scale'    Either 'on' to scale each column of X by its standard deviation
%                before fitting, or 'off' (the default) to omit scaling.
%
%   Example:
%      load hald 
%      stepwisefit(ingredients,heat,'penter',.08)

%   Copyright 1993-2015 The MathWorks, Inc.


narginchk(2,Inf);

okargs   = {'inmodel' 'penter' 'premove' 'display' 'maxiter' 'keep' 'scale'};
defaults = {[]        []       []        'on'      Inf       []    'off'};
[in,penter,premove,dodisp,maxiter,keep,scale] = ...
                internal.stats.parseArgs(okargs,defaults,varargin{:});

P = size(allx,2);
if isempty(in)
   in = false(1,P);
elseif islogical(in)
   if length(in)~=P
      error(message('stats:stepwisefit:BadInModel'));
   end
else
   if any(~ismember(in,1:P))
      error(message('stats:stepwisefit:BadInModel'));
   end
   in = ismember((1:P),in);
end

if isempty(keep)
   keep = false(size(in));
elseif islogical(keep)
   if length(keep)~=P
      error(message('stats:stepwisefit:BadKeep'));
   end
else
   if any(~ismember(keep,1:P))
      error(message('stats:stepwisefit:BadKeep'));
   end
   keep = ismember((1:P),keep);
end

% Get p-to-enter and p-to-remove defaults
if isempty(penter) && isempty(premove)
   penter = 0.05;
   premove = 0.10;
elseif isempty(penter)
   penter = min(premove,0.05);
elseif isempty(premove)
   premove = max(penter,0.10);
end
if numel(penter)~=1 || ~isnumeric(penter) || penter<=0 || penter>=1
   error(message('stats:stepwisefit:BadPEnter'));
end
if numel(premove)~=1 || ~isnumeric(premove) || premove<=0 || premove>=1
   error(message('stats:stepwisefit:BadPRemove'));
end
if penter>premove
   error(message('stats:stepwisefit:BadPEnterRemove'));
end

% Check input dimensions
if size(y,2)~=1
   error(message('stats:stepwisefit:InvalidData'));
end
if size(y,1)~=size(allx,1)
   error(message('stats:stepwisefit:InputSizeMismatch'));
end
if ~isreal(allx) || ~isreal(y)
    error(message('stats:stepwisefit:ComplexData'));
end
   

% Remove NaN rows, if any
if any(any(isnan(allx))) || any(any(isnan(y)))
   [badin,wasnan,allx,y] = statremovenan(allx,y);
   if (badin>0)
      error(message('stats:stepwisefit:InputSizeMismatch'))
   end
else
   wasnan = false(size(y));
end

% Determine if results are to be presented based on standardize X columns,
% but standardize internally in all cases
doscale = isequal(scale,'on');
sx = std(allx,0,1);
sx(sx==0) = 1;
allx = allx./sx(ones(size(allx,1),1),:);
   
% If requested, display information about the starting state
if isequal(dodisp,'on')
   if ~any(in)
      coltext = 'none';
   else
      coltext = sprintf('%d ',find(in));
   end      
   fprintf(getString(message('stats:stepwisefit:cltext_InitCols',sprintf('%s\n',coltext))));
end

% Set up variables that describe the step history
if nargout>=7
   Bhist = zeros(length(sx),0);
   rmse = [];
   df0 = [];
   inmat = false(0,length(in));
end

% Start iteratively moving terms in and out
jstep = 0;
while(true)
   % Perform current fit
   [B,SE,PVAL,stats] = stepcalc(allx,y,in,id,t,approach);
   if ~doscale
       % Undo scaling if this was not requested
       B = B./sx';
       SE = SE./sx';
       stats.B = B;
       stats.SE = SE;
       stats.xr = stats.xr .* repmat(sx(:,~in),size(stats.xr,1),1);
       stats.covb = stats.covb ./ (sx'*sx);
   end

   % Select next step
   [nextstep,pinout] = stepnext(in,PVAL,B,penter,premove,keep);

   % Remember what happened in this step
   if nargout>=7 && jstep>0
      Bhist(in,jstep) = B(in);
      rmse(jstep) = stats.rmse;
      df0(jstep) = stats.df0;
      inmat(jstep,:) = in;
   end


   if (jstep>=maxiter), break; end
   jstep = jstep + 1;

   % Report the action for this step
   if nextstep==0
      break;
   elseif isequal(dodisp,'on')
      if in(nextstep)
         disp(getString(message('stats:stepwisefit:cltext_NextStepRemoved', ...
             sprintf('%d',jstep), sprintf('%d',nextstep), sprintf('%g',pinout))));
      else
         disp(getString(message('stats:stepwisefit:cltext_NextStepAdded', ...
             sprintf('%d',jstep), sprintf('%d',nextstep), sprintf('%g',pinout))));
      end
   end
   in(nextstep) = ~in(nextstep);
end

if isequal(dodisp,'on')
   if ~any(in)
      coltext = 'none';
   else
      coltext = sprintf('%d ',find(in));
   end      
   fprintf(getString(message('stats:stepwisefit:cltext_FinalCols', ...
       sprintf('%s\n',coltext))));
   inout = {getString(message('stats:stepwisefit:cltext_Out')); ...
       getString(message('stats:stepwisefit:cltext_In'))};
    disp([{getString(message('stats:stepwisefit:cltext_Coeff')) ...
           getString(message('stats:stepwisefit:cltext_StdErr')) ...
           getString(message('stats:stepwisefit:cltext_Status')) ...
           getString(message('stats:stepwisefit:cltext_Pval'))}; ...
           num2cell(B) num2cell(SE) inout(in+1) num2cell(PVAL)]);
end

% Remember which rows were removed
if nargout>=5
   stats.wasnan = wasnan;
end

% Return history of steps taken
if nargout>=7
   history.B = Bhist;
   history.rmse = rmse;
   history.df0 = df0;
   history.in = inmat;
end

% -----------------------------------------
function [swap,p] = stepnext(in,PVAL,B,penter,premove,keep)
%STEPNEXT Figure out next step

swap = 0;
p = NaN;

% Look for terms out that should be in
termsout = find(~in & ~keep);
if ~isempty(termsout)
   [pmin,kmin] = min(PVAL(termsout));
   if pmin<penter
      swap = termsout(kmin(1));
      p = pmin;
   end
end

% Otherwise look for terms in that should be out
if swap==0
   termsin = find(in & ~keep);
   if ~isempty(termsin)
      badterms = termsin(isnan(PVAL(termsin)));
      if ~isempty(badterms)
          % Apparently we have a perfect fit but it is also overdetermined.
          % Terms with NaN coefficients may as well be removed.
          swap = isnan(B(badterms));
          if any(swap)
              swap = badterms(swap);
              swap = swap(1);
          else              
              % If there are many terms contributing to a perfect fit, we
              % may as well remove the term that contributes the least.
              % For convenience we'll pick the one with the smallest coeff.
              [~,swap] = min(abs(B(badterms)));
              swap = badterms(swap);
          end
          p = NaN;
      else
          [pmax,kmax] = max(PVAL(termsin));
          if pmax>premove
             swap = termsin(kmax(1));
             p = pmax;
          end
      end
   end
end


% -----------------------------------------
function [B,SE,PVAL,stats] = stepcalc(allx,y,in,id,t,approach)
%STEPCALC Perform fit and other calculations as part of stepwise regression



N = length(y);
P = length(in);
X = [ones(N,1) allx(:,in)];
nin = sum(in)+1;
tol = max(N,P+1)*eps(class(allx));
x = allx(:,~in);
sumxsq = sum(x.^2,1);

% Rrank = size(X,2);

% Compute b and its standard error
[Q,R,perm] = qr(X,0);
if isempty(R)
    Rrank = 0;
else
    Rrank = sum(abs(diag(R)) > tol*abs(R(1)));
end
if Rrank < nin
    R = R(1:Rrank,1:Rrank);
    Q = Q(:,1:Rrank);
    perm = perm(1:Rrank);
end

% Compute the LS coefficients, filling in zeros in elements corresponding
% to rows of X that were thrown out.
b = zeros(nin,1);
Qb = Q'*y;
Qb(abs(Qb) < tol*max(abs(Qb))) = 0;
b(perm) = R \ Qb;

if strcmp('markov',approach)
    [~,~,results] = qls(id,y,t,X,1,approach);
else
    [~,~,results] = gee(id,y,t,X,1,approach);
end
b = [results.robust{3:end,2}]';

r = y - X*b;
dfe = size(X,1)-Rrank;
df0 = Rrank - 1;
SStotal = norm(y-mean(y))^2;
SSresid = norm(r)^2;
perfectyfit = (dfe==0) || (SSresid<tol*SStotal);
if perfectyfit
    SSresid = 0;
    r = 0*r;
end
rmse = sqrt(safedivide(SSresid,dfe));
Rinv = R\eye(size(R));
se = zeros(nin,1);
covb = zeros(nin);
covb(perm,perm) = rmse^2 * (Rinv*Rinv');
se = sqrt(diag(covb));

% Compute separate added-variable coeffs and their standard errors
xr = x - Q*(Q'*x);  % remove effect of "in" predictors on "out" predictors
yr = r;             % remove effect of "in" predictors on response

xx = sum(xr.^2,1);

perfectxfit = (xx<=tol*sumxsq);
if any(perfectxfit)      % to coef==0 for columns dependent in "in" cols
    xr(:,perfectxfit) = 0;  
    xx(perfectxfit) = 1;
end
b2 = safedivide(yr'*xr, xx);
r2 = repmat(yr,1,sum(~in)) - xr .* repmat(b2,N,1);
df2 = max(0,dfe - 1);
s2 = safedivide(sqrt(safedivide(sum(r2.^2,1),df2)), sqrt(xx));

% Combine in/out coefficients and standard errors
B = zeros(P,1);
B(in) = b(2:end);
B(~in) = b2';
SE = zeros(P,1);
SE(in) = se(2:end);
SE(~in) = s2';
COVB = zeros(P,P);
COVB(in,in) = covb(2:end,2:end);

% Get P-to-enter or P-to-remove for each term
PVAL = zeros(P,1);
tstat = zeros(P,1);
if any(in)
   tval = safedivide(B(in),SE(in));
   ptemp = 2*tcdf(-abs(tval),dfe);
   PVAL(in) = ptemp;
   tstat(in) = tval;
end
if any(~in)
   if dfe>1
      tval = safedivide(B(~in),SE(~in));
      ptemp = 2*tcdf(-abs(tval),dfe-1);
   else
      tval = NaN;
      ptemp = NaN;
   end
   PVAL(~in) = ptemp;
   tstat(~in) = tval;
end

% Compute some summary statistics
MSexplained = safedivide(SStotal-SSresid, df0);
fstat = safedivide(MSexplained, rmse^2);
pval = fpval(fstat,df0,dfe);

% Return summary statistics as a single structure
stats.source = 'stepwisefit';
stats.dfe = dfe;
stats.df0 = df0;
stats.SStotal = SStotal;
stats.SSresid = SSresid;
stats.fstat = fstat;
stats.pval = pval;
stats.rmse = rmse;
stats.xr = xr;
stats.yr = yr;
stats.B = B;
stats.SE = SE;
stats.TSTAT = tstat;
stats.PVAL = PVAL;
stats.covb = COVB;
stats.intercept = b(1);

% --------------------------------------
function quotient = safedivide(num,denom)
t = (denom==0);
if ~any(t) || isempty(num)
    quotient = num ./ denom;
else
    if isscalar(num) && ~isscalar(denom)
        num = repmat(num,size(denom));
    elseif isscalar(denom) && ~isscalar(num)
        denom = repmat(denom,size(num));
        t = (denom==0);
    end
    quotient(~t) = num(~t) ./ denom(~t);
    quotient(t) = Inf * sign(num(t));
end
