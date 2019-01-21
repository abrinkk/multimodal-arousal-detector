function SD = postprocess2(ar_p,w_p)
%POSTPROCESS2 combines arousal and wake predictions as a single measure.
%   SD = POSTPROCESS2(ar_p,w_p) computes the union of ar_p and w_p,
%   whereafter predictions closer than 10 seconds are merged. 
%
%   Author: Andreas Brink-Kjaer.
%   Date: 17-Jun-2018
%
%   Input:  ar_p, arousal predictions
%           w_p, wake predictions
%   Output: SD, combined arousal measure

% Union of ar_p and w_p
SD = ar_p | w_p;
% Connect < 10 s
bw_non_sd = bwlabel(~SD);
for i = 2:max(bw_non_sd)-1
    if sum(bw_non_sd == i) < 10
        SD(bw_non_sd == i) = 1;
    end
end
end