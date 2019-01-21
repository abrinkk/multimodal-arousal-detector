function w2 = w_postprocess(w,T)
%W_POSTPROCESS processes wake probability to get predictions.
%   w2 = W_POSTPROCESS(w,T) performs postprocessing of wake
%   probability by threshold and duration rules to get wake predictions.
%
%   Author: Andreas Brink-Kjaer.
%   Date: 17-Jun-2018
%
%   Input:  w, wake probability in 1 second bins
%           T, threshold value [0, 1]
%   Output: w2, wake predictions in 1 second bins

% Threshold
w = w > T;
% Connect < 15 s
bw_non_w = bwlabel(~w);
for i = 2:max(bw_non_w)-1
    if sum(bw_non_w == i) < 15
        w(bw_non_w == i) = 1;
    end
end
% Remove < 15 seconds wake
bw_w = bwlabel(w);
for i = 1:max(bw_w)
    if sum(bw_w == i) < 15
        w(bw_w == i) = 0;
    end
end
w2 = w;
end