function e = anc_rls_m(sig,ref1,ref2)
% function e = anc_rls_m(sig,ref1,ref2)
% Adaptive noise canceller based on recursive least squares algorithm as
% developed/expounded on in P. He paper, "Removal of ocular artifracts from
% electro-encephalogram by adaptive filtering" - Medical and Biological
% Engineering and Computing 2004, Vol. 42
%
% e = cleaned output signal
% sig = input signal contaminated by noise from reference signals ref1 and ref2
%
%
% written by Hyatt Moore IV, February 28, 2012
%

if(nargin==3)
    %make sure we are dealing with row vectors
    [r,c] =size(sig);
    if(r<c)
        sig = sig.';
        ref1 = ref1.';
        ref2 = ref2.';
    end
    
    %(a) set initial values
    lambda = 0.995; %forgetting factor; He. found values between 0.995 and 1.0 to produce similarly good result
    lambda_inv = 1/lambda;
    M = 4; %number of weights for the filter; why 4?  because I like the number 4, and He. found that after 3 things are fine.
    sigma = 0.01;
    Rinv = eye(2*M)/sigma;   %Rinv is the [R(n-1)]^-1
    H = zeros(M*2,1);  %contains weights of both filters (hv and hh = horizontal and vertical)
    
    e = zeros(size(sig));
    
    for n=M:numel(sig)
        r = [ref1(n-M+1:n);ref2(n-M+1:n)];
        K = (Rinv*r)/(lambda+r'*Rinv*r);
        e_n = sig(n)-r'*H;
        H = H+K*e_n;
        Rinv = lambda_inv*Rinv-lambda_inv*K*r'*Rinv;
        e(n) = sig(n)-r'*H;
    end
    
elseif(nargin==2)
    
    %make sure we are dealing with row vectors
    [r,c] =size(sig);
    
    if(r<c)
        sig = sig.';
        ref1 = ref1.';
    end
    
    assert(all(size(sig) == size(ref1)));
    
    %(a) set initial values
    lambda = 0.995; %forgetting factor; He. found values between 0.995 and 1.0 to produce similarly good result
    lambda_inv = 1/lambda;
    M = 4; %number of weights for the filter; why 4?  because I like the number 4, and He. found that after 3 things are fine.
    sigma = 0.01;
    Rinv = eye(M)/sigma;   %Rinv is the [R(n-1)]^-1
    H = zeros(M,1);  %contains weights of both filters (hv and hh = horizontal and vertical)
    
    e = zeros(size(sig));
    
    for n=M:numel(sig)
        r = ref1(n-M+1:n);
        K = (Rinv*r)/(lambda+r'*Rinv*r);
        e_n = sig(n)-r'*H;
        H = H+K*e_n;
        Rinv = lambda_inv*Rinv-lambda_inv*K*r'*Rinv;
        e(n) = sig(n)-r'*H;
    end
end