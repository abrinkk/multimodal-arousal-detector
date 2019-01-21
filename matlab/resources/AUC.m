function [PR_x, PR_y, auPR] = AUC(Precision, Recall)
%AUC calculates the area under a curve.
%   [PR_x, PR_y, auPR] = AUC(Precision, Recall) calculates the
%   area under the precision-recall curve. 
%
%   Author: Caspar Aleksander Bang Jespersen
%   Modified by Andreas Brink-Kjaer
%   Date: 17-Jun-2018
%
%   Input:  Precision, precision scores for different thresholds.
%           Recall, recall scores for different thresholds.
%   Output: PR_x, sorted recall scores
%           PR_y, sorted precision scores
%           auPR, area under precision recall curve.

% Remove NANs
valid = ~(isnan(Recall) | isnan(Precision));
Recall = Recall(valid);
Precision = Precision(valid);

% Add boundary conditions
x_ = [0; Recall; 1];
y_ = [1; Precision; 0];

% Order x axis
[~,order] = sort(x_);
x_ = x_(order);
y_ = y_(order);

% Remove duplicates
xu = unique(x_);
xh = hist(x_,xu);
duplicates = find(xh>1);
for i = 1:length(duplicates)
    idx = find(x_==xu(duplicates(i)));
    ynew = mean(y_(idx));
    y_(idx(1)) = ynew;
    x_(idx(2:end)) = nan;
end
valid = ~isnan(x_);
x_ = x_(valid);
y_ = y_(valid);

% Interpolate
PR_x = min(x_):0.01:1;
% PR_y = zeros(size(PR_x));
% for i = 1:length(PR_x)
%     if any(PR_x(i) == x_)
%         PR_y(i) = x_(PR_x(i) == x_);
%     else
%         idx_1 = find(x_ < PR_x(i),1);
%         idx_2 = find(x_ > PR_x(i),1);
%         TPx = round(PR_x(i)*(TP1 + FN1) - TP1);
%         PR_y(i) = (TP1 + TPx)/(TP1 + TPx + FP1 + (FP2 - FP1)/(TP2 - TP1)*TPx);
%     end
% end
%         
PR_y = interp1(x_,y_,PR_x,'linear');
auPR = mean(PR_y);