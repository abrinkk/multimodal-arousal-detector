function [ idx_seconds, idx_samples ] = time2ind(T, fs )
%TIME2IND converts HH:MM:SS.FF format to seconds from record start.
%   [ idx_seconds, idx_samples ] = TIME2IND(T, fs ) inputs event table from
%   WSC2 and converts the time in HH:MM:SS.FF format in T.Time to seconds.
%
%   Author: Caspar Aleksander Bang Jespersen
%   Modified by Andreas Brink-Kjaer
%   Date: 17-Jun-2018
%
%   Input:  T, table from WSC2 annotation file
%           fs, sampling frequency
%   Output: idx_seconds, event times in seconds
%           idd_samples, event times in samples

if ~exist('fs','var')
    fs = 128;
end
timenum = datenum(T.Time);

% Events to start
start_first = find(ismember(T.Event,'Start Recording'));
start_resume = find(ismember(T.Event,'Recording Resumed'));
stop_pause = find(ismember(T.Event,'Paused'));
start_events = sort([start_first; start_resume]);
stop_events = stop_pause;

% Assert they add up
assert(length(start_events) == length(stop_events));

% Collect sessions
s = cell(length(start_events), 2);
s(:,1) = T.Time(start_events);
s(:,2) = T.Time(stop_events);

% Convert to datenums and make up for change in day (restrict
% that they must be monotonically increasing)
s = cellfun(@datenum, s)';
for i = 2:numel(s)
    if s(i) < s(i-1)
        s(i) = s(i) + 1;
    end
end
s = s';

% Adjust input times for date change
timenum(timenum < s(1)) = timenum(timenum < s(1)) + 1;

% Find time instance of event
cds = cumsum(s(:,2)-s(:,1));
dds = [0; cds(1:end-1)];
try
    idx_s = arrayfun(@(num) find((num >= s(:,1)) & (num <= s(:,2))), timenum,'Un',0);
    % For scorings inside pauses
    idx_inpause = cellfun(@(c) isempty(c),idx_s);
    if any(idx_inpause)
        idx_s(idx_inpause) = idx_s(find(idx_inpause)+1);
        timenum(idx_inpause) = s([idx_s{idx_inpause}],1);
    end
    idx_s = cell2mat(idx_s);
    % Set times
    idx_day = dds(idx_s) + (timenum - s(idx_s,1));
    idx_seconds = (idx_day) * (24*60*60);
    idx_samples = 1 + round((idx_seconds) * (fs));
catch
    idx_s = 0;
    idx_day = 0;
    idx_seconds = 0;
    idx_samples = 0;
end
end
