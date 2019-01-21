function y = filterEKGArtefact(sig_data, ref_data)
%PREPROCESS.FILTEREKGARTEFACT performs RLS adaptive filtering.
%   y = PREPROCESS.FILTEREKGARTEFACT(sig_data, ref_data) performs adaptive 
%   noise cancelling based on recursive least squares algorithm. 
%
%   Author: Andreas Brink-Kjaer.
%   Date: 17-Jun-2018
%
%   Input:  sig_data, input contaminated data (EEG, EOG, or EMG)
%           ref_data, noise signal (ECG)
%   Output: y, cleaned data

% Use shortest signal
least = min([length(sig_data) length(ref_data)]);
% Adaptive filter
e = filter.anc_rls_m(sig_data(1:least)', ref_data(1:least)');
% Remove possible NAN for flat segments
e(isnan(e)) = sig_data(isnan(e));
% Insert filtered signal
y = zeros(size(sig_data));
y(1:length(e)) = e;
y = y(:)';
end