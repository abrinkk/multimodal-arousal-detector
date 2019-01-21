%PREPROCESSDATA process inputs edfs and labels
%   PREPROCESSDATA inputs all edfs and labels for used databases and
%   filters, resamples and normalizes data. The data is now saved in a new
%   .txt format.
%
%   Author: Andreas Brink-Kjaer.
%   Date: 17-Jun-2018
%

clear all; close all;
des_fs = 128;
dirIndex = paths;
trainRatio = 0.7;
valRatio = 0.1;
testRatio = 0.2;
data_paths;

%% Iteratre through all data
Overwrite = 0;

%% MrOS
parfor i = 1:length(f_mros_edf)
    fprintf('Processsing MrOS %.0f/%.0f\n',i,length(f_mros_edf));
    try
        if i/length(f_mros_edf) < trainRatio
            folder = 'train';
        elseif i/length(f_mros_edf) < trainRatio + valRatio
            folder = 'val';
        else
            folder = 'test'; 
        end
        if ~Overwrite && exist(filepath(dirIndex.Data,folder,[f_mros_edf{i}(1:end-4) '.txt']),'file')
            continue;
        end
        % Load and preprocess data
        [hdr,data] = LoadEDF(filepath(p_mros_edf,f_mros_edf{i}));
        [hdr,data] = preprocess.resampleData(data,hdr,des_fs);
        data = preprocess.filterData(data,hdr,1);
        
        % Load labels
        [~,ar_seq] = LoadAR(filepath(p_mros_lab,f_mros_lab{i}),size(data,2)/des_fs);
        W = LoadWake(filepath(p_mros_lab,f_mros_lab{i}),size(data,2)/des_fs);
        
        % Save
        exportData(data,ar_seq,W,hdr,filepath(dirIndex.Data,folder,[f_mros_edf{i}(1:end-4) '.txt']));
    catch me
        disp(me.message);
    end
end

%% CFS
parfor i = 1:length(f_cfs_edf)
    fprintf('Processsing CFS %.0f/%.0f\n',i,length(f_cfs_edf));
    try
        if i/length(f_cfs_edf) < trainRatio
            folder = 'train';
        elseif i/length(f_cfs_edf) < trainRatio + valRatio
            folder = 'val';
        else
            folder = 'test';
        end
        if ~Overwrite && exist(filepath(dirIndex.Data,folder,[f_cfs_edf{i}(1:end-4) '.txt']),'file')
            continue;
        end
        % Load and preprocess data
        [hdr,data] = LoadEDF(filepath(p_cfs_edf,f_cfs_edf{i}));
        [hdr,data] = preprocess.resampleData(data,hdr,des_fs);
        data = preprocess.filterData(data,hdr,1);
        
        % Load labels
        [~,ar_seq] = LoadAR(filepath(p_cfs_lab,f_cfs_lab{i}),size(data,2)/des_fs);
        W = LoadWake(filepath(p_cfs_lab,f_cfs_lab{i}),size(data,2)/des_fs);
        
        % Save
        exportData(data,ar_seq,W,hdr,filepath(dirIndex.Data,folder,[f_cfs_edf{i}(1:end-4) '.txt']));
    catch me
        disp(me.message);
    end
end

%% WSC2
parfor i = 1:length(f_wsc2_edf)
    fprintf('Processsing WSC2 %.0f/%.0f\n',i,length(f_wsc2_edf));
    try
        folder = 'other';
        if ~Overwrite && exist(filepath(dirIndex.Data,folder,[f_wsc2_edf{i}(1:end-4) '.txt']),'file')
            continue;
        end
        % Load and preprocess data
        [hdr,data] = LoadEDF(filepath(p_wsc2_edf,f_wsc2_edf{i}),'wsc2');
        [hdr,data] = preprocess.resampleData(data,hdr,des_fs);
        data = preprocess.filterData(data,hdr,1);
        
        % Load labels
        [~,ar_seq] = LoadAR(filepath(p_wsc2_lab,f_wsc2_lab{i}),size(data,2)/des_fs,'wsc2');
        W = LoadWake(filepath(p_wsc2_lab,f_wsc2_lab{i}),size(data,2)/des_fs,'wsc2');
        
        % Save
        exportData(data,ar_seq,W,hdr,filepath(dirIndex.Data,folder,[f_wsc2_edf{i}(1:end-4) '.txt']));
    catch me
        disp(me.message);
    end
end

%% SSC
parfor i = 1:length(f_ssc_edf)
    fprintf('Processsing SSC %.0f/%.0f\n',i,length(f_ssc_edf));
    try
        folder = 'other';
        if ~Overwrite && exist(filepath(dirIndex.Data,folder,[f_ssc_edf{i}(1:end-4) '.txt']),'file')
            continue;
        end
        % Load and preprocess data
        [hdr,data] = LoadEDF(filepath(p_ssc_edf,f_ssc_edf{i}),'ssc');
        [hdr,data] = preprocess.resampleData(data,hdr,des_fs);
        data = preprocess.filterData(data,hdr,1);
        
        % Load labels
        [~,ar_seq] = LoadAR(filepath(p_ssc_lab,f_ssc_lab{i}),size(data,2)/des_fs,'ssc');
        W = LoadWake(filepath(p_ssc_lab,f_ssc_lab{i}),size(data,2)/des_fs,'ssc');
        
        % Save
        exportData(data,ar_seq,W,hdr,filepath(dirIndex.Data,folder,[f_ssc_edf{i}(1:end-4) '.txt']));
    catch me
        disp(me.message);
    end
end

%% WSC
parfor i = 1:length(f_wsc_edf)
    fprintf('Processsing WSC %.0f/%.0f\n',i,length(f_wsc_edf));
    try
        folder = 'other';
        if ~Overwrite && exist(filepath(dirIndex.Data,folder,[f_wsc_edf{i}(1:end-4) '.txt']),'file')
            continue;
        end
        % Load and preprocess data
        [hdr,data] = LoadEDF(filepath(p_wsc_edf,f_wsc_edf{i}),'wsc');
        [hdr,data] = preprocess.resampleData(data,hdr,des_fs);
        data = preprocess.filterData(data,hdr,1);
        
        % Load labels
        [~,ar_seq] = LoadAR(filepath(p_wsc_lab,f_wsc_lab{i}),size(data,2)/des_fs,'wsc');
        W = LoadWake(filepath(p_wsc_lab,f_wsc_lab{i}),size(data,2)/des_fs,'wsc');
        
        % Save
        exportData(data,ar_seq,W,hdr,filepath(dirIndex.Data,folder,[f_wsc_edf{i}(1:end-4) '.txt']));
    catch me
        disp(me.message);
    end
end