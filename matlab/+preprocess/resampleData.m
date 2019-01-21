function [hdr,dataR] = resampleData(data,hdr,des_fs)
%PREPROCESS.RESAMPLEDATA resamples data to a desired samplings frequency.
%   [hdr,dataR] = PREPROCESS.RESAMPLEDATA(data,hdr,des_fs) resamples the
%   input data to a desired sampling frequency using an anti-alias filter.
%
%   Author: Andreas Brink-Kjaer.
%   Date: 17-Jun-2018
%
%   Input:  data, EEG, EOG, EMG and ECG signals extracted with LoadEDF.
%           hdr, .edf header extracted with LoadEDF.
%           des_fs, desired sampling frequency.
%   Output: hdr, .edf header with altered sampling frequencies
%           dataR, resampled data

% Data length
max_fs = max(hdr.fs);
L = ceil(size(data,2)/max_fs*des_fs);
dataR = nan(size(data,1),L);
% Iterate over channels
for i = 1:size(data,1)
    % End of signal in data structure
    dataEnd = ceil(hdr.fs(i)/max_fs*size(data,2));
    dataClip = data(i,1:dataEnd);
    % Check up/down sample factor
    [p,q] = rat(des_fs/hdr.fs(i));
    if hdr.fs(i) ~= des_fs
        % Resample
        dataResample = resample(dataClip,p,q);
        dataR(i,1:round(dataEnd*p/q)) = dataResample;
    else
        dataR(i,1:dataEnd) = dataClip;
    end
end
% Change header sampling frequency
hdr.fs = ones(size(hdr.fs))*des_fs;
end
