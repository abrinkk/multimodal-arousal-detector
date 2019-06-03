function [dirIndex, profile] = paths()
%PATHS add required paths and file directories.
%   [dirIndex, profile] = PATHS() add paths to functions used.
%
%   Author: Caspar Aleksander Bang Jespersen.
%   Modified by Andreas Brink-Kjaer.
%   Date: 17-Jun-2018
%
%   Input:  None
%   Output: dirIndex, a struct containing directories of commonly used
%           folders
%           profile, the used computing profile, e.g. local or sherlock
%           server.

persistent pthDirIndex pthProfile;

if isempty(pthDirIndex) || isempty(pthProfile)
    
    pthDirIndex = struct;
    
    if ~exist('profile','var')
        if exist('../profile','file')
            pthProfile = cellstr(fileread('../profile'));
        else
            pthProfile = 'abk';
        end
    end
    if ismember('new',pthProfile)
        pathtofile = mfilename('fullpath');
        pthDirIndex.base = pathtofile(1:end-12);
        pthDirIndex.matlabCode = [pthDirIndex.base 'matlab\'];
        pthDirIndex.Data = [pthDirIndex.matlabCode 'data\'];
        pthDirIndex.architecture = 'i64';
        pthDirIndex.python = '';
    elseif ismember('abk',pthProfile)
        pthDirIndex.base = 'C:\Users\andre\Desktop\ardetector2018\';
        pthDirIndex.matlabCode = 'C:\Users\andre\Desktop\ardetector2018\matlab\';
        pthDirIndex.Data = [pthDirIndex.matlabCode 'data\'];
        pthDirIndex.architecture = 'i64';
        pthDirIndex.python = '';
    elseif ismember('abk2',pthProfile)
        pthDirIndex.base = 'C:\Users\andre\OneDrive\Skrivebord\ardetector2018\';
        pthDirIndex.matlabCode = 'C:\Users\andre\OneDrive\Skrivebord\ardetector2018\matlab\';
        pthDirIndex.Data = 'D:\ardetector\data\';
        pthDirIndex.architecture = 'i64';
        pthDirIndex.python = '';
    elseif ismember('sherlock',pthProfile)
        pthDirIndex.base = '/scratch/users/abk26/ardetector-master/';
        pthDirIndex.matlabCode = '/scratch/users/abk26/ardetector-master/matlab/';
        pthDirIndex.Data = '/scratch/users/abk26/ardetector-master/matlab/data/';
        pthDirIndex.architecture = 'a64';
        pthDirIndex.python = 'python';
    elseif ismember('mignotserver',pthProfile)
        pthDirIndex.base = '/data1/home/abk26/ardetector2018/';
        pthDirIndex.matlabCode = '/data1/home/abk26/ardetector2018/matlab/';
        pthDirIndex.Data = '/data1/home/abk26/ardetector2018/matlab/data/';
        pthDirIndex.architecture = 'a64';
        pthDirIndex.python = 'python';
    else
        error('Profile not found.');
    end
    % Add paths
    addpath([pthDirIndex.matlabCode 'resources_matlab']);
    addpath(genpath([pthDirIndex.matlabCode 'resources_matlab/edfread']));
    addpath(genpath([pthDirIndex.matlabCode 'resources_matlab/export_fig']));
    addpath(genpath([pthDirIndex.matlabCode 'resources_matlab/append_pdfs']));
    addpath(genpath([pthDirIndex.matlabCode 'resources_matlab/GEEQBOX']));
    
end

dirIndex = pthDirIndex;
profile = pthProfile;