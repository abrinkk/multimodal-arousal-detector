function [hdr,data] = LoadEDF(p_file,ftype)
%LOADEDF reads edf files and selects desired channels.
%   [hdr,data] = LOADEDF(p_file,ftype) loads the edf file specified in
%   p_file. Depending on the .edf file source, a set of channels are
%   extracted and referenced to match system specifications.
%
%   Author: Andreas Brink-Kjaer.
%   Date: 17-Jun-2018
%
%   Input:  p_file, .edf file location
%           ftype, string of data source
%   Output: hdr, edf file header
%           data, EEG, EOG, EMG and ECG channels

% Determine data type
if ~exist('ftype','var')
    if contains(p_file,'mros')
        ftype = 'mros';
    elseif contains(p_file,'cfs')
        ftype = 'cfs';
    end
end
% Load data
[hdr,data] = edfread(p_file);

switch ftype
    case 'mros'
        % Get signal derivations C3/A2, ROC/A1, LOC/A2, LChin/RChin, ECGL/ECGR.
        A1 = contains(hdr.label,'A1');
        A2 = contains(hdr.label,'A2');
        C3 = contains(hdr.label,'C3');
        LOC = contains(hdr.label,'LOC');
        ROC = contains(hdr.label,'ROC');
        EMG1 = contains(hdr.label,'RChin');
        EMG2 = contains(hdr.label,'LChin');
        ECG2 = contains(hdr.label,'ECGL');
        ECG1 = contains(hdr.label,'ECGR');
        EEG = data(C3,:) -  data(A2,:);
        EOGR = data(ROC,:) -  data(A1,:);
        EOGL = data(LOC,:) -  data(A2,:);
        EMG = data(EMG2,:) -  data(EMG1,:);
        ECG = data(ECG2,:) -  data(ECG1,:);
        data = [EEG; EOGR; EOGL; EMG; ECG];
        % Index header
        idx = [find(C3) find(ROC) find(LOC) find(EMG1) find(ECG2)];
        hdr.label = {'EEGC','EOGR','EOGL','EMG','ECG'};
        hdr.ns = 5;
        hdr.transducer = hdr.transducer(idx);
        hdr.units = hdr.units(idx);
        hdr.physicalMin = hdr.physicalMin(idx);
        hdr.physicalMax = hdr.physicalMax(idx);
        hdr.digitalMin = hdr.digitalMin(idx);
        hdr.digitalMax = hdr.digitalMax(idx);
        hdr.prefilter = hdr.prefilter(idx);
        hdr.samples = hdr.samples(idx);
        hdr.fs = hdr.samples./hdr.duration;
        for i = 1:length(hdr.units)
            if strcmp(hdr.units{i},'mV')
                data(i,:) = data(i,:)*1000;
                hdr.units{i} = 'uV';
            end
        end
    case 'cfs'
        % Get signal derivations C3/M2, ROC/M1, LOC/M2, EMG1/EMG2, ECG2/ECG1.
        M1 = contains(hdr.label,'M1');
        M2 = contains(hdr.label,'M2');
        C3 = contains(hdr.label,'C3');
        LOC = contains(hdr.label,'LOC');
        ROC = contains(hdr.label,'ROC');
        EMG1 = contains(hdr.label,'EMG1');
        EMG2 = contains(hdr.label,'EMG2');
        ECG2 = contains(hdr.label,'ECG2');
        ECG1 = contains(hdr.label,'ECG1');
        EEG = data(C3,:) -  data(M2,:);
        EOGR = data(ROC,:) -  data(M1,:);
        EOGL = data(LOC,:) -  data(M2,:);
        EMG = data(EMG2,:) -  data(EMG1,:);
        ECG = data(ECG2,:) -  data(ECG1,:);
        data = [EEG; EOGR; EOGL; EMG; ECG];
        % Index header
        idx = [find(C3) find(ROC) find(LOC) find(EMG1) find(ECG2)];
        data = data(:,1:max(hdr.samples(idx))/max(hdr.samples)*size(data,2));
        hdr.label = {'EEGC','EOGR','EOGL','EMG','ECG'};
        hdr.ns = 5;
        hdr.transducer = hdr.transducer(idx);
        hdr.units = hdr.units(idx);
        hdr.physicalMin = hdr.physicalMin(idx);
        hdr.physicalMax = hdr.physicalMax(idx);
        hdr.digitalMin = hdr.digitalMin(idx);
        hdr.digitalMax = hdr.digitalMax(idx);
        hdr.prefilter = hdr.prefilter(idx);
        hdr.samples = hdr.samples(idx);
        hdr.fs = hdr.samples./hdr.duration;
        for i = 1:length(hdr.units)
            if strcmp(hdr.units{i},'mV')
                data(i,:) = data(i,:)*1000;
                hdr.units{i} = 'uV';
            end
        end
    case {'wsc2','ssc','wsc'}
        % Get signal derivations
        Channel_EEG_Cen = {'C3M2' 'C3M1' 'C4M2' 'C4M1' 'C4AVG' ...
            'C3AVG','C3A2','C4A1','C3A1','C4A2','C3x','C4x'};
        Channel_REOG = {'REOGM1','ROCM1','ROCA1','ROCA2','REOGx'};
        Channel_LEOG = {'LEOGM2','LOCM1','LOCA2','LOCA1','LEOGx'};
        Channel_EMG = {'Chin1Chin2' 'Chin1Chin3' 'Chin3Chin2' 'ChinEMG'};
        Channel_ECG = {'ECG' 'EKG1AVG' 'EKG1EKG2' 'LLEG1EKG2' 'EKG'};
        EEG_idx = find(contains(hdr.label,Channel_EEG_Cen),1,'first');
        EOGR_idx = find(contains(hdr.label,Channel_REOG),1,'first');
        EOGL_idx = find(contains(hdr.label,Channel_LEOG),1,'first');
        EMG_idx = find(contains(hdr.label,Channel_EMG),1,'first');
        ECG_idx = find(contains(hdr.label,Channel_ECG),1,'first');
        EEG = data(EEG_idx,:);
        EOGR = data(EOGR_idx,:);
        EOGL = data(EOGL_idx,:);
        EMG = data(EMG_idx,:);
        ECG = data(ECG_idx,:);
        data = [EEG; EOGR; EOGL; EMG; ECG];
        idx = [EEG_idx EOGR_idx EOGL_idx EMG_idx ECG_idx];
        hdr.label = {'EEGC','EOGR','EOGL','EMG','ECG'};
        hdr.ns = 5;
        hdr.transducer = hdr.transducer(idx);
        hdr.units = hdr.units(idx);
        hdr.physicalMin = hdr.physicalMin(idx);
        hdr.physicalMax = hdr.physicalMax(idx);
        hdr.digitalMin = hdr.digitalMin(idx);
        hdr.digitalMax = hdr.digitalMax(idx);
        hdr.prefilter = hdr.prefilter(idx);
        hdr.samples = hdr.samples(idx);
        hdr.fs = hdr.samples./hdr.duration;
        for i = 1:length(hdr.units)
            if strcmp(hdr.units{i},'mV')
                data(i,:) = data(i,:)*1000;
                hdr.units{i} = 'uV';
            end
        end
end
end
