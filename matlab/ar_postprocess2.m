function ar2 = ar_postprocess2(ar,w,Tar,Tw)
% Threshold
ar = ar > Tar;
w = w > Tw;
% Connect ar < 10 s
bw_non_ar = bwlabel(~ar);
for i = 2:max(bw_non_ar)-1
    if sum(bw_non_ar == i) < 10
        ar(bw_non_ar == i) = 1;
    end
end
% Connect w < 15 s
bw_non_w = bwlabel(~w);
for i = 2:max(bw_non_ar)-1
    if sum(bw_non_ar == i) < 15
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