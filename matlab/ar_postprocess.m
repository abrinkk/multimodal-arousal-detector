function ar2 = ar_postprocess(ar,T,remove3)
%AR_POSTPROCESS processes arousal probability to get predictions.
%   ar2 = AR_POSTPROCESS(ar,T) performs postprocessing of arousal
%   probability by threshold and duration rules to get arousal predictions.
%   This postprocessing also removes short arousals in match AASM rules,
%   thereby enabling a more accurate comparison to the gold standard.
%
%   Author: Andreas Brink-Kjaer.
%   Date: 17-Jun-2018
%
%   Input:  ar, arousal probability in 1 second bins
%           T, threshold value [0, 1]
%           remove3, boolean to remove < 3 second arousals (default = 1).
%   Output: ar2, arousal predictions in 1 second bins

if ~exist('remove3','var')
    remove3 = 1;
end

% Threshold
ar = ar > T;
% Connect < 10 s
bw_non_ar = bwlabel(~ar);
for i = 2:max(bw_non_ar)-1
    if sum(bw_non_ar == i) < 10
        ar(bw_non_ar == i) = 1;
    end
end
% Remove < 3 s seconds arousals
if remove3
    bw_ar = bwlabel(ar);
    for i = 1:max(bw_ar)
        if sum(bw_ar == i) < 3
            ar(bw_ar == i) = 0;
        end
    end
end
ar2 = ar;
end