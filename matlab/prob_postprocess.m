function [ar,w] = prob_postprocess(ar_prob,w_prob,Tar,Tw)
% Threshold
ar = ar_prob > Tar;
w = w_prob > Tw;
% Connect < 10 s
bw_non_ar = bwlabel(~ar);
for i = 2:max(bw_non_ar)-1
    if sum(bw_non_ar == i) < 10
        ar(bw_non_ar == i) = 1;
    end
end
% Remove < 3 s seconds arousals
bw_ar = bwlabel(ar);
for i = 1:max(bw_ar)
    if sum(bw_ar == i) < 3
        ar(bw_ar == i) = 0;
    end
end
ar2 = ar;
end