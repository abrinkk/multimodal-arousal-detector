function [ar_prob,w_prob,ar_prob2,w_prob2] = getPred(fname,T1,T2,L,bs)
%GETPRED reads model predictions and performs postprocessing.
%   [ar_prob,w_prob] = GETPRED(fname,T1,T2,L) loads model predictions and
%   performs postprocessing of arousal and wake probabilities.
%
%   Author: Andreas Brink-Kjaer.
%   Date: 17-Jun-2018
%   
%   Input:  fname, filename
%           T1, arousal postprocessing threshold
%           T2, wake postprocessing threshold
%           L, prediction length
%   Output: ar_prob, arousal predictions
%           w_prob, wake predictions
%           ar_prob2, arousal prediction with batch edge error
%           w_prob2, wake predictions with batch edge error

if ~exist('T1','var')
    T1 = 0;
end
if ~exist('T2','var')
    T2 = 0;
end
if ~exist('L','var')
    L = 0;
end


% Read prediction files
FID = fopen(fname,'r');
C = textscan(FID,'%s');
C = C{1};
ar_prob = textscan(C{1},'%.2f','Delimiter',',');
ar_prob = ar_prob{1};
w_prob = textscan(C{2},'%.2f','Delimiter',',');
w_prob = w_prob{1};
ar_prob2 = textscan(C{3},'%.2f','Delimiter',',');
ar_prob2 = ar_prob2{1};
w_prob2 = textscan(C{4},'%.2f','Delimiter',',');
w_prob2 = w_prob2{1};
fclose(FID);

if ~exist('bs','var')
    bs = length(ar_prob) - length(ar_prob2);
end

% Add empty batch
ar_prob2 = [zeros(bs/2,1); ar_prob2; zeros(bs/2,1)];
w_prob2 = [zeros(bs/2,1); w_prob2; zeros(bs/2,1)];

% Temporary
ar_temp = ar_prob;
w_temp = w_prob;

% Combine predictions
nb = length(ar_prob)/bs;
idx2 = sort(unique(reshape(bs*(1:(nb-1))' + (-bs/4:(bs/4-1)),1,[])));
ar_prob(idx2) = ar_prob2(idx2);
w_prob(idx2) = w_prob2(idx2);

% Save raw batching
ar_prob2 = ar_temp;
w_prob2 = w_temp;

% Postprocessing if specified
if T1 > 0
    ar_prob = ar_postprocess(ar_prob,T1);
end
if T2 > 0
    w_prob = w_postprocess(w_prob,T2);
end

% Add length
if L > 0 && L ~= length(ar_prob)
    ar_prob = [ar_prob; zeros(L-length(ar_prob),1)];
    w_prob = [w_prob; ones(L-length(w_prob),1)];
    ar_prob2 = [ar_prob2; zeros(L-length(ar_prob2),1)];
    w_prob2 = [w_prob2; ones(L-length(w_prob2),1)];
end
end