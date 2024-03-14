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
    else
        ftype = 'default';
    end
end

% Load data
if strcmp(ftype, 'default')
    [data, hdr] = lab_read_edf2(p_file,[],0);
    hdr = hdr.hdr;
    hdr.label = strtrim(cellstr(hdr.channelname)');
    hdr.label = erase(hdr.label,[" ","-"]);
elseif strcmp(ftype,'kassel')
    [data,hdr] = lab_read_edf(p_file);
    hdr = hdr.hdr;
    hdr.label = strtrim(cellstr(hdr.channelname)');
    hdr.label = erase(hdr.label,[" ","-"]);
elseif strcmp(ftype, 'stages') || contains(ftype, 'dcsm') || strcmp(ftype, 'mayo_ond') || strcmp(ftype, 'iceberg')
    [data, hdr] = lab_read_edf2(p_file,[],0);
    hdr = hdr.hdr;
    hdr.label = strtrim(cellstr(hdr.channelname)');
    hdr.label = erase(hdr.label,[" ","-"]);
elseif strcmp(ftype, 'dreem')
    h5_info = h5info(p_file);
    new_channels = {h5_info.Datasets.Name}';
elseif strcmp(ftype, 'jedi')
    h5_info = h5info(p_file);
else
    [hdr,data] = edfread(p_file);
end

switch ftype
    case 'default'
        % Should consider all cases in english with standard naming
        [data, hdr] = default_edf_loader(data, hdr);
    case 'mayo_ond'
        [data, hdr] = default_edf_loader(data, hdr);
    case 'iceberg'
        % Should consider all cases in english with standard naming
        channel_struct = struct('REOG', {get_channel_alias(construct_channel_name({'R', 'EOG', 'E2', 'ROC', 'REOG', 'EOGR'}, {'D', 'R', 'EOG', '2', 'A1', 'A2', 'M1', 'M2', 'A11', 'A12', 'A21', 'A22', 'M11', 'M12', 'M21', 'M22', 'x', 'avg'}, {'', '-'}))}, ...
                                'LEOG', {get_channel_alias(construct_channel_name({'L', 'EOG', 'E1', 'LOC', 'LEOG', 'EOGL'}, {'G', 'L', 'EOG', '2', 'A1', 'A2', 'M1', 'M2', 'A11', 'A12', 'A21', 'A22', 'M11', 'M12', 'M21', 'M22', 'x', 'avg'}, {'', '-'}))});
        [data, hdr] = default_edf_loader(data, hdr, channel_struct);
        
    case {'dcsm_rbd', 'dcsm_ptsd'}
        % Add alternate names
        alt_names = @(x) unique([x, strrep(strrep(strrep(strrep(strrep(strrep(x,' ',''),'-',''),'#',''),'(',''),')',''),'_','')]);
        % Get signal derivations
        % Referenced
        Channel_EEG_Cen = alt_names({'C3-M2','C3M2',...
            'EEG C3-A2','EEG C3-A22','C3AVG','C3A2','C3A1','C3x','C3:A2','c3-a2','EEG C3-A2'});
        Channel_REOG = alt_names({'R-EOG','EOG2','E2M2','EOG ROC-A2','EOG ROC-A22'...
            'REOGM1','ROCM1','ROCA1','ROCA2','REOGx','EOGr:A2','eogr-a1','EOGH-A1','EEG EOGH-A1','EEG EOGH-A1-Ref','E2'});
        Channel_LEOG = alt_names({'L-EOG','EOG1','E1M2','EOG LOC-A1','EOG LOC-A2','EOG LOC-A22','EOG LOC-A11'...
            'E1M2','LEOGM2','LOCM1','LOCA2','LOCA1','LEOGx','EOGl:A2','eogl-a2','EOGL-A2','EEG EOGV-A2','EOGV-A2','EEG EOGV-A2-Ref','E1'});
        Channel_EMG = alt_names({'EMG Chin','EMG Aux12','CHIN','Chin',...
            'Chin1Chin2' 'Chin1Chin3' 'Chin3Chin2' 'ChinEMG','EMG','CHIN','Chin','chin','EMG Chin','EMG CHIN','EMG','Chin-Gnd'});
        Channel_ECG = alt_names({'ECG IIHF','EKG',...
            'ECG' 'EKG1AVG' 'EKG1EKG2' 'LLEG1EKG2' 'EKG','ECG2','ECG','ECG EKG','ecg','ECG II', 'ECG 2','EKG','EKG-Gnd','EMG EKG'});
        Channel_C3 = alt_names({'C3','EEG C3-Ref','C3-Ref'});
        Channel_ROC = alt_names({'E2 (REOG)','E2(REOG)','E2','ROC', 'EOGr','EOGH','EEG EOGH-Ref','EOGH-Ref'});
        Channel_LOC = alt_names({'E1 (LEOG)','E1(LEOG)','E1','LOC', 'EOGV','EEG EOGV-Ref','EOGV-Ref','EOGl', 'EOGr'});
        Channel_A2 = alt_names({'M2','A2','A2-Ref','EEG A2-Ref'});
        Channel_A1 = alt_names({'M1','A1','A1-Ref','EEG A1-Ref'});
        Channel_EMG2 = alt_names({'EMG #2','EMG Aux2','EMG Chin2','EMG2','CHIN2','Chin 2','Chin2','EMG2','EMG+','ChinR'});
        Channel_EMG1 = alt_names({'EMG #1','EMG Aux1','EMG1','CHIN1','Chin 1','Chin1','EMG-','ChinL'});
        Channel_ECG2 = alt_names({'ECG 2','ECG I2','ECG II','ECG II2','ECG2','EKG #2','EKG2'});
        Channel_ECG1 = alt_names({'ECG 1','ECG I','ECG1','EKG #1','EKG1'});
        
        % Referneced
        EEG_idx = find(ismember(hdr.label,Channel_EEG_Cen),1,'first');
        EOGR_idx = find(ismember(hdr.label,Channel_REOG),1,'first');
        EOGL_idx = find(ismember(hdr.label,Channel_LEOG),1,'first');
        EMG_idx = find(ismember(hdr.label,Channel_EMG),1,'first');
        ECG_idx = find(ismember(hdr.label,Channel_ECG),1,'first');
        % Unreferenced
        C3_idx = find(ismember(hdr.label,Channel_C3),1,'first');
        ROC_idx = find(ismember(hdr.label,Channel_ROC),1,'first');
        LOC_idx = find(ismember(hdr.label,Channel_LOC),1,'first');
        A1_idx = find(ismember(hdr.label,Channel_A1),1,'first');
        A2_idx = find(ismember(hdr.label,Channel_A2),1,'first');
        EMG1_idx = find(ismember(hdr.label,Channel_EMG1),1,'first');
        EMG2_idx = find(ismember(hdr.label,Channel_EMG2),1,'first');
        ECG1_idx = find(ismember(hdr.label,Channel_ECG1),1,'first');
        ECG2_idx = find(ismember(hdr.label,Channel_ECG2),1,'first');
        % Check if exist
        if ~isempty(EEG_idx)
            EEG = data{EEG_idx};
        else
            EEG = data{C3_idx} - data{A2_idx};
            EEG_idx = C3_idx;
        end
        if ~isempty(EOGR_idx)
            EOGR = data{EOGR_idx};
        else
            EOGR = data{ROC_idx} - data{A2_idx};
            EOGR_idx = ROC_idx;
        end
        if ~isempty(EOGL_idx)
            EOGL = data{EOGL_idx};
        else
            EOGL = data{LOC_idx} - data{A1_idx};
            EOGL_idx = LOC_idx;
        end
        if ~isempty(EMG_idx)
            EMG = data{EMG_idx};
        else
            EMG = data{EMG2_idx} - data{EMG1_idx};
            EMG_idx = EMG2_idx;
        end
        if ~isempty(ECG_idx)
            ECG = data{ECG_idx};
            no_ECG = false;
        elseif ~isempty(ECG2_idx)
            ECG = data{ECG2_idx} - data{ECG1_idx};
            ECG_idx = ECG2_idx;
            no_ECG = false;
        else
            no_ECG = true;
            ECG = 0;
            ECG_idx = EMG_idx;
        end
        max_dur = max([length(EEG), length(EOGR), length(EOGL), length(EMG), length(ECG)]);
        EEG = [EEG, zeros(1,max_dur - length(EEG))];
        EOGR = [EOGR, zeros(1,max_dur - length(EOGR))];
        EOGL = [EOGL, zeros(1,max_dur - length(EOGL))];
        EMG = [EMG, zeros(1,max_dur - length(EMG))];
        ECG = [ECG, zeros(1,max_dur - length(ECG))];
        data = [EEG; EOGR; EOGL; EMG; ECG];
        idx = [EEG_idx EOGR_idx EOGL_idx EMG_idx ECG_idx];
        hdr.label = {'EEGC','EOGR','EOGL','EMG','ECG'};
        hdr.ns = 5;
        hdr.transducer = hdr.transducer(idx);
%         hdr.units = hdr.units(idx);
%         hdr.physicalMin = hdr.physicalMin(idx);
%         hdr.physicalMax = hdr.physicalMax(idx);
%         hdr.digitalMin = hdr.digitalMin(idx);
%         hdr.digitalMax = hdr.digitalMax(idx);
%         hdr.prefilter = hdr.prefilter(idx);
        hdr.fs = hdr.numbersperrecord(idx) / hdr.duration;
%         hdr.samples = hdr.samples(idx);
%         hdr.fs = hdr.samples./hdr.duration;
%         for i = 1:length(hdr.units)
%             if strcmp(hdr.units{i},'mV')
%                 data(i,:) = data(i,:)*1000;
%                 hdr.units{i} = 'uV';
%             end
%         end
        hdr.no_ECG = no_ECG;
    case 'umhan'
        % Add alternate names
        alt_names = @(x) unique([x, strrep(strrep(strrep(strrep(strrep(strrep(x,' ',''),'-',''),'#',''),'(',''),')',''),'_','')]);
        % Get signal derivations
        % Referenced
        Channel_EEG_Cen = alt_names({'C3-M2','C3M2',...
            'EEG C3-A2','EEG C3-A22','C3AVG','C3A2','C3A1','C3x'});
        Channel_REOG = alt_names({'R-EOG','EOG2','E2M1','E2M2','EOG ROC-A2','EOG ROC-A22'...
            'REOGM1','ROCM1','ROCA1','ROCA2','REOGx','EOGrA2'});
        Channel_LEOG = alt_names({'L-EOG','EOG1','E1M2','EOG LOC-A1','EOG LOC-A2','EOG LOC-A22','EOG LOC-A11'...
            'E1M2','LEOGM2','LOCM1','LOCA2','LOCA1','LEOGx','EOGlA2'});
        Channel_EMG = alt_names({'EMG Chin','EMG Aux12','CHIN','Chin',...
            'Chin1Chin2' 'Chin1Chin3' 'Chin3Chin2' 'ChinEMG', 'EMG'});
        Channel_ECG = alt_names({'ECG IIHF','EKG',...
            'ECG' 'EKG1AVG' 'EKG1EKG2' 'LLEG1EKG2' 'EKG'});
        Channel_C3 = alt_names({'C3'});
        Channel_ROC = alt_names({'E2 (REOG)','E2(REOG)','E2','ROC'});
        Channel_LOC = alt_names({'E1 (LEOG)','E1(LEOG)','E1','LOC'});
        Channel_A2 = alt_names({'M2','A2'});
        Channel_A1 = alt_names({'M1','A1'});
        Channel_EMG2 = alt_names({'EMG #2','EMG Aux2','EMG Chin2','EMG2','CHIN2','Chin 2','Chin2'});
        Channel_EMG1 = alt_names({'EMG #1','EMG Aux1','EMG1','CHIN1','Chin 1','Chin1'});
        Channel_ECG2 = alt_names({'ECG 2','ECG I2','ECG II','ECG II2','ECG2','EKG #2','EKG2'});
        Channel_ECG1 = alt_names({'ECG 1','ECG I','ECG1','EKG #1','EKG1'});
        
        % Referneced
        EEG_idx = find(ismember(hdr.label,Channel_EEG_Cen),1,'first');
        EOGR_idx = find(ismember(hdr.label,Channel_REOG),1,'first');
        EOGL_idx = find(ismember(hdr.label,Channel_LEOG),1,'first');
        EMG_idx = find(ismember(hdr.label,Channel_EMG),1,'first');
        ECG_idx = find(ismember(hdr.label,Channel_ECG),1,'first');
        
        % Unreferenced
        C3_idx = find(ismember(hdr.label,Channel_C3),1,'first');
        ROC_idx = find(ismember(hdr.label,Channel_ROC),1,'first');
        LOC_idx = find(ismember(hdr.label,Channel_LOC),1,'first');
        A1_idx = find(ismember(hdr.label,Channel_A1),1,'first');
        A2_idx = find(ismember(hdr.label,Channel_A2),1,'first');
        EMG1_idx = find(ismember(hdr.label,Channel_EMG1),1,'first');
        EMG2_idx = find(ismember(hdr.label,Channel_EMG2),1,'first');
        ECG1_idx = find(ismember(hdr.label,Channel_ECG1),1,'first');
        ECG2_idx = find(ismember(hdr.label,Channel_ECG2),1,'first');
        % Check if exist
        if ~isempty(EEG_idx)
            EEG = data{EEG_idx};
        else
            EEG = data{C3_idx} - data{A2_idx};
            EEG_idx = C3_idx;
        end
        if ~isempty(EOGR_idx)
            EOGR = data{EOGR_idx};
        else
            EOGR = data{ROC_idx} - data{A2_idx};
            EOGR_idx = ROC_idx;
        end
        if ~isempty(EOGL_idx)
            EOGL = data{EOGL_idx};
        else
            EOGL = data{LOC_idx} - data{A1_idx};
            EOGL_idx = LOC_idx;
        end
        if ~isempty(EMG_idx)
            EMG = data{EMG_idx};
        elseif ~isempty(EMG2_idx)
            EMG = data{EMG2_idx} - data{EMG1_idx};
            EMG_idx = EMG2_idx;
        else
            EMG = [1, -1, zeros(1, length(EEG) - 2)];
            EMG_idx = EEG_idx;
        end
        if ~isempty(ECG_idx)
            ECG = data{ECG_idx};
            no_ECG = false;
        elseif ~isempty(ECG2_idx)
            ECG = data{ECG2_idx} - data{ECG1_idx};
            ECG_idx = ECG2_idx;
            no_ECG = false;
        else
            no_ECG = true;
            ECG = 0;
            ECG_idx = EMG_idx;
        end
        max_dur = max([length(EEG), length(EOGR), length(EOGL), length(EMG), length(ECG)]);
        EEG = [EEG, zeros(1,max_dur - length(EEG))];
        EOGR = [EOGR, zeros(1,max_dur - length(EOGR))];
        EOGL = [EOGL, zeros(1,max_dur - length(EOGL))];
        EMG = [EMG, zeros(1,max_dur - length(EMG))];
        ECG = [ECG, zeros(1,max_dur - length(ECG))];
        data = [EEG; EOGR; EOGL; EMG; ECG];
        idx = [EEG_idx EOGR_idx EOGL_idx EMG_idx ECG_idx];
        hdr.label = {'EEGC','EOGR','EOGL','EMG','ECG'};
        hdr.ns = 5;
        hdr.transducer = hdr.transducer(idx);
%         hdr.units = hdr.units(idx);
%         hdr.physicalMin = hdr.physicalMin(idx);
%         hdr.physicalMax = hdr.physicalMax(idx);
%         hdr.digitalMin = hdr.digitalMin(idx);
%         hdr.digitalMax = hdr.digitalMax(idx);
%         hdr.prefilter = hdr.prefilter(idx);
        hdr.fs = hdr.numbersperrecord(idx);
%         hdr.samples = hdr.samples(idx);
%         hdr.fs = hdr.samples./hdr.duration;
%         for i = 1:length(hdr.units)
%             if strcmp(hdr.units{i},'mV')
%                 data(i,:) = data(i,:)*1000;
%                 hdr.units{i} = 'uV';
%             end
%         end
        hdr.no_ECG = no_ECG;
    case 'multi'
        alt_names = @(x) unique([x, strrep(strrep(strrep(strrep(strrep(strrep(x,' ',''),'-',''),'#',''),'(',''),')',''),'_','')]);
        % Get signal derivations
        Channel_EEG_Cen = alt_names({'C3-A2'});
        Channel_REOG = alt_names({'ROC-A1','ROC-A2'});
        Channel_LEOG = alt_names({'LOC-A2'});
        Channel_EMG = alt_names({'EMG1-EMG2'});
        Channel_ECG = alt_names({'ECG1-ECG2'});
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
    case 'jedi'
        new_channels = {h5_info.Groups(strcmp({h5_info.Groups.Name}, '/signals')).Datasets.Name};
        fs = double(h5_info.Attributes(strcmp({h5_info.Attributes.Name}, 'sampling_rate')).Value);
        units = h5_info.Attributes(strcmp({h5_info.Attributes.Name}, 'unit_voltage')).Value;
        
        % Get signal derivations
        Channel_EEG_Cen = {'c3a2','c3m2','c3-m2','c3-a2'};
        Channel_REOG = {'e2','eog2a1','e2-m1','e2-a1'};
        Channel_LEOG = {'e1','eog1a2','e1-m2','e1-a2'};
        Channel_EMG = {'chin1-chin2','chin','chin2'};
        Channel_ECG = {'ecg','ecg2','ecg_2'};
        EEG_n = Channel_EEG_Cen{find(ismember(Channel_EEG_Cen,new_channels),1,'first')};
        EOGR_n = Channel_REOG{find(ismember(Channel_REOG,new_channels),1,'first')};
        EOGL_n = Channel_LEOG{find(ismember(Channel_LEOG,new_channels),1,'first')};
        EMG_n = Channel_EMG{find(ismember(Channel_EMG,new_channels),1,'first')};
        ECG_n = Channel_ECG{find(ismember(Channel_ECG,new_channels),1,'first')};
        EEG_idx = find(ismember(new_channels, EEG_n));
        EOGR_idx = find(ismember(new_channels, EOGR_n));
        EOGL_idx = find(ismember(new_channels, EOGL_n));
        EMG_idx = find(ismember(new_channels, EMG_n));
        ECG_idx = find(ismember(new_channels, ECG_n));
        EEG = h5read(p_file,['/signals/' EEG_n]);
        EOGR = h5read(p_file,['/signals/' EOGR_n]);
        EOGL = h5read(p_file,['/signals/' EOGL_n]);
        EMG = h5read(p_file,['/signals/' EMG_n]);
        ECG = h5read(p_file,['/signals/' ECG_n]);
        % TODO: ZERO PAD FIRST
        max_dur = max([length(EEG), length(EOGR), length(EOGL), length(EMG), length(ECG)]);
        EEG = [EEG, zeros(1,max_dur - length(EEG))];
        EOGR = [EOGR, zeros(1,max_dur - length(EOGR))];
        EOGL = [EOGL, zeros(1,max_dur - length(EOGL))];
        EMG = [EMG, zeros(1,max_dur - length(EMG))];
        ECG = [ECG, zeros(1,max_dur - length(ECG))];
        data = double([EEG; EOGR; EOGL; EMG; ECG]); 
        idx = [EEG_idx EOGR_idx EOGL_idx EMG_idx ECG_idx];
        hdr = struct;
        hdr.label = {'EEGC','EOGR','EOGL','EMG','ECG'};
        hdr.ns = 5;
        
        if iscell(units)
            units = units{0};
        end
        hdr.units = arrayfun(@(x) units, 1:hdr.ns, 'Un', 0);
        hdr.fs = repmat(fs, 1, 5);
        for i = 1:length(hdr.units)
            if strcmp(hdr.units{i},'V')
                data(i,:) = data(i,:)*1000*1000;
                hdr.units{i} = 'uV';
            end
            if strcmp(hdr.units{i},'mV')
                data(i,:) = data(i,:)*1000;
                hdr.units{i} = 'uV';
            end
        end
    case 'dreem'
        get_chan_attributes = @(x, attr) h5_info.Datasets(x).Attributes(strcmp({h5_info.Datasets(x).Attributes.Name},attr)).Value;
        % Get signal derivations
        Channel_EEG_Cen = {'c3a2','c3m2'};
        Channel_REOG = {'e2','eog2a1'};
        Channel_LEOG = {'e1','eog1a2'};
        Channel_EMG = {'chinchin2','chin','chin2'};
        Channel_ECG = {'ecg','ecg2','ecg_2'};
        EEG_n = Channel_EEG_Cen{find(ismember(Channel_EEG_Cen,new_channels),1,'first')};
        EOGR_n = Channel_REOG{find(ismember(Channel_REOG,new_channels),1,'first')};
        EOGL_n = Channel_LEOG{find(ismember(Channel_LEOG,new_channels),1,'first')};
        EMG_n = Channel_EMG{find(ismember(Channel_EMG,new_channels),1,'first')};
        ECG_n = Channel_ECG{find(ismember(Channel_ECG,new_channels),1,'first')};
        EEG_idx = find(ismember(new_channels, EEG_n));
        EOGR_idx = find(ismember(new_channels, EOGR_n));
        EOGL_idx = find(ismember(new_channels, EOGL_n));
        EMG_idx = find(ismember(new_channels, EMG_n));
        ECG_idx = find(ismember(new_channels, ECG_n));
        EEG = h5read(p_file,['/' EEG_n])';
        EOGR = h5read(p_file,['/' EOGR_n])';
        EOGL = h5read(p_file,['/' EOGL_n])';
        EMG = h5read(p_file,['/' EMG_n])';
        ECG = h5read(p_file,['/' ECG_n])';
        % TODO: ZERO PAD FIRST
        max_dur = max([length(EEG), length(EOGR), length(EOGL), length(EMG), length(ECG)]);
        EEG = [EEG, zeros(1,max_dur - length(EEG))];
        EOGR = [EOGR, zeros(1,max_dur - length(EOGR))];
        EOGL = [EOGL, zeros(1,max_dur - length(EOGL))];
        EMG = [EMG, zeros(1,max_dur - length(EMG))];
        ECG = [ECG, zeros(1,max_dur - length(ECG))];
        data = double([EEG; EOGR; EOGL; EMG; ECG]); 
        idx = [EEG_idx EOGR_idx EOGL_idx EMG_idx ECG_idx];
        hdr = struct;
        hdr.label = {'EEGC','EOGR','EOGL','EMG','ECG'};
        hdr.ns = 5;
        
        hdr.units = arrayfun(@(x) char(get_chan_attributes(x, 'physical_dimension')), idx, 'Un',0);
        hdr.fs = arrayfun(@(x) double(get_chan_attributes(x, 'frequency')), idx);
        for i = 1:length(hdr.units)
            if strcmp(hdr.units{i},'mV')
                data(i,:) = data(i,:)*1000;
                hdr.units{i} = 'uV';
            end
        end
    case 'stages'
        % Add alternate names
        alt_names = @(x) unique([x, strrep(strrep(strrep(strrep(strrep(strrep(x,' ',''),'-',''),'#',''),'(',''),')',''),'_','')]);
        % Get signal derivations
        % Referenced
        Channel_EEG_Cen = alt_names({'C3-M2','C3M2',...
            'EEG C3-A2','EEG C3-A22','C3AVG','C3A2','C3A1','C3x'});
        Channel_REOG = alt_names({'R-EOG','EOG2','E2M2','EOG ROC-A2','EOG ROC-A22'...
            'REOGM1','ROCM1','ROCA1','ROCA2','REOGx'});
        Channel_LEOG = alt_names({'L-EOG','EOG1','E1M2','EOG LOC-A1','EOG LOC-A2','EOG LOC-A22','EOG LOC-A11'...
            'E1M2','LEOGM2','LOCM1','LOCA2','LOCA1','LEOGx'});
        Channel_EMG = alt_names({'EMG Chin','EMG Aux12','CHIN','Chin',...
            'Chin1Chin2' 'Chin1Chin3' 'Chin3Chin2' 'ChinEMG'});
        Channel_ECG = alt_names({'ECG IIHF','EKG',...
            'ECG' 'EKG1AVG' 'EKG1EKG2' 'LLEG1EKG2' 'EKG'});
        Channel_C3 = alt_names({'C3'});
        Channel_ROC = alt_names({'E2 (REOG)','E2(REOG)','E2','ROC'});
        Channel_LOC = alt_names({'E1 (LEOG)','E1(LEOG)','E1','LOC'});
        Channel_A2 = alt_names({'M2','A2'});
        Channel_A1 = alt_names({'M1','A1'});
        Channel_EMG2 = alt_names({'EMG #2','EMG Aux2','EMG Chin2','EMG2','CHIN2','Chin 2','Chin2'});
        Channel_EMG1 = alt_names({'EMG #1','EMG Aux1','EMG1','CHIN1','Chin 1','Chin1'});
        Channel_ECG2 = alt_names({'ECG 2','ECG I2','ECG II','ECG II2','ECG2','EKG #2','EKG2'});
        Channel_ECG1 = alt_names({'ECG 1','ECG I','ECG1','EKG #1','EKG1'});
        if contains(p_file,'STLK')
            Channel_REOG = alt_names({'R-EOG','EOG2','E2M2','EOG ROC-A2','EOG ROC-A22'...
                'REOGM1','ROCM1','ROCA1','ROCA2','REOGx','E2'});
            Channel_LEOG = alt_names({'L-EOG','EOG1','E1M2','EOG LOC-A1','EOG LOC-A2','EOG LOC-A22','EOG LOC-A11'...
                'E1M2','LEOGM2','LOCM1','LOCA2','LOCA1','LEOGx','E1'});
        end
        
        % Referneced
        EEG_idx = find(ismember(hdr.label,Channel_EEG_Cen),1,'first');
        EOGR_idx = find(ismember(hdr.label,Channel_REOG),1,'first');
        EOGL_idx = find(ismember(hdr.label,Channel_LEOG),1,'first');
        EMG_idx = find(ismember(hdr.label,Channel_EMG),1,'first');
        ECG_idx = find(ismember(hdr.label,Channel_ECG),1,'first');
        % Unreferenced
        C3_idx = find(ismember(hdr.label,Channel_C3),1,'first');
        ROC_idx = find(ismember(hdr.label,Channel_ROC),1,'first');
        LOC_idx = find(ismember(hdr.label,Channel_LOC),1,'first');
        A1_idx = find(ismember(hdr.label,Channel_A1),1,'first');
        A2_idx = find(ismember(hdr.label,Channel_A2),1,'first');
        EMG1_idx = find(ismember(hdr.label,Channel_EMG1),1,'first');
        EMG2_idx = find(ismember(hdr.label,Channel_EMG2),1,'first');
        ECG1_idx = find(ismember(hdr.label,Channel_ECG1),1,'first');
        ECG2_idx = find(ismember(hdr.label,Channel_ECG2),1,'first');
        % Check if exist
        if ~isempty(EEG_idx)
            EEG = data{EEG_idx};
        else
            EEG = data{C3_idx} - data{A2_idx};
            EEG_idx = C3_idx;
        end
        if ~isempty(EOGR_idx)
            EOGR = data{EOGR_idx};
        else
            EOGR = data{ROC_idx} - data{A2_idx};
            EOGR_idx = ROC_idx;
        end
        if ~isempty(EOGL_idx)
            EOGL = data{EOGL_idx};
        else
            EOGL = data{LOC_idx} - data{A1_idx};
            EOGL_idx = LOC_idx;
        end
        if ~isempty(EMG_idx)
            EMG = data{EMG_idx};
        else
            EMG = data{EMG2_idx} - data{EMG1_idx};
            EMG_idx = EMG2_idx;
        end
        if ~isempty(ECG_idx)
            ECG = data{ECG_idx};
            no_ECG = false;
        elseif ~isempty(ECG2_idx)
            ECG = data{ECG2_idx} - data{ECG1_idx};
            ECG_idx = ECG2_idx;
            no_ECG = false;
        else
            no_ECG = true;
            ECG = 0;
            ECG_idx = EMG_idx;
        end
        max_dur = max([length(EEG), length(EOGR), length(EOGL), length(EMG), length(ECG)]);
        EEG = [EEG, zeros(1,max_dur - length(EEG))];
        EOGR = [EOGR, zeros(1,max_dur - length(EOGR))];
        EOGL = [EOGL, zeros(1,max_dur - length(EOGL))];
        EMG = [EMG, zeros(1,max_dur - length(EMG))];
        ECG = [ECG, zeros(1,max_dur - length(ECG))];
        data = [EEG; EOGR; EOGL; EMG; ECG];
        idx = [EEG_idx EOGR_idx EOGL_idx EMG_idx ECG_idx];
        hdr.label = {'EEGC','EOGR','EOGL','EMG','ECG'};
        hdr.ns = 5;
        hdr.transducer = hdr.transducer(idx);
%         hdr.units = hdr.units(idx);
%         hdr.physicalMin = hdr.physicalMin(idx);
%         hdr.physicalMax = hdr.physicalMax(idx);
%         hdr.digitalMin = hdr.digitalMin(idx);
%         hdr.digitalMax = hdr.digitalMax(idx);
%         hdr.prefilter = hdr.prefilter(idx);
        hdr.fs = hdr.numbersperrecord(idx);
%         hdr.samples = hdr.samples(idx);
%         hdr.fs = hdr.samples./hdr.duration;
%         for i = 1:length(hdr.units)
%             if strcmp(hdr.units{i},'mV')
%                 data(i,:) = data(i,:)*1000;
%                 hdr.units{i} = 'uV';
%             end
%         end
        hdr.no_ECG = no_ECG;
    case 'mesa'
        % Get signal derivations C3/A2, ROC/A1, LOC/A2, LChin/RChin, ECGL/ECGR.
        C3 = contains(hdr.label,'EEG3');
        ROC = contains(hdr.label,'EOGR');
        LOC = contains(hdr.label,'EOGL');
        EMG1 = contains(hdr.label,'EMG');
        ECG2 = contains(hdr.label,'EKG');
        EEG = data(C3,:);
        EOGR = data(ROC,:);
        EOGL = data(LOC,:);
        EMG = data(EMG1,:);
        ECG = data(ECG2,:);
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
    case 'shhs'
        % Get signal derivations C3/A2, ROC/A1, LOC/A2, LChin/RChin, ECGL/ECGR.
        C3 = contains(hdr.label,'EEGsec');
        ROC = contains(hdr.label,'EOGR');
        LOC = contains(hdr.label,'EOGL');
        EMG1 = contains(hdr.label,'EMG');
        ECG2 = contains(hdr.label,'ECG');
        EEG = data(C3,:);
        EOGR = data(ROC,:);
        EOGL = data(LOC,:);
        EMG = data(EMG1,:);
        ECG = data(ECG2,:);
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
    case 'kassel'
        % if unreferenced
        if all(~contains(hdr.label,'EEGC3A2'))
            % Get signal derivations C3/A2, ROC/A1, LOC/A2, LChin/RChin, ECGL/ECGR.
            A1 = contains(hdr.label,'A1');
            A2 = contains(hdr.label,'A2');
            C3 = strcmp(hdr.label,'C3');
            LOC = contains(hdr.label,'LOC');
            ROC = contains(hdr.label,'ROC');
            EMG1 = contains(hdr.label,'CHIN1');
            EMG2 = contains(hdr.label,'CHIN2');
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
%             hdr.units = hdr.units(idx);
            hdr.physicalMin = hdr.physmin(idx);
            hdr.physicalMax = hdr.physmax(idx);
            hdr.digitalMin = hdr.digimin(idx);
            hdr.digitalMax = hdr.digimax(idx);
            hdr.prefilter = hdr.prefilt(idx);
            hdr.samples = hdr.numbersperrecord(idx);
            hdr.fs = hdr.samples./hdr.duration;
%             for i = 1:length(hdr.units)
%                 if strcmp(hdr.units{i},'mV')
%                     data(i,:) = data(i,:)*1000;
%                     hdr.units{i} = 'uV';
%                 end
%             end
        else
            EEG_idx = contains(hdr.label,'EEGC3A2');
            EOGR_idx = contains(hdr.label,'EEGROCA1');
            EOGL_idx = contains(hdr.label,'EEGLOCA2');
            EMG_idx = contains(hdr.label,'EMGCHIN1CHIN2');
            ECG_idx = contains(hdr.label,'ECGECGLECGR');
            EEG = data(EEG_idx,:);
            EOGR = data(EOGR_idx,:);
            EOGL = data(EOGL_idx,:);
            EMG = data(EMG_idx,:);
            ECG = data(ECG_idx,:);
            data = [EEG; EOGR; EOGL; EMG; ECG];
            idx = [find(EEG_idx) find(EOGR_idx) find(EOGL_idx) find(EMG_idx) find(ECG_idx)];
            hdr.label = {'EEGC','EOGR','EOGL','EMG','ECG'};
            hdr.ns = 5;
            hdr.transducer = hdr.transducer(idx);
%             hdr.units = hdr.units(idx);
            hdr.physicalMin = hdr.physmin(idx);
            hdr.physicalMax = hdr.physmax(idx);
            hdr.digitalMin = hdr.digimin(idx);
            hdr.digitalMax = hdr.digimax(idx);
            hdr.prefilter = hdr.prefilt(idx);
            hdr.samples = hdr.numbersperrecord(idx);
            hdr.fs = hdr.samples./hdr.duration;
%             for i = 1:length(hdr.units)
%                 if strcmp(hdr.units{i},'mV')
%                     data(i,:) = data(i,:)*1000;
%                     hdr.units{i} = 'uV';
%                 end
%             end
        end
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
    case 'chat'
        % Get signal derivations C3/M2, ROC/M1, LOC/M2, EMG1/EMG2, ECG2/ECG1.
        M1 = contains(hdr.label,'M1','IgnoreCase',1);
        M2 = contains(hdr.label,'M2','IgnoreCase',1);
        C3 = contains(hdr.label,'C3','IgnoreCase',1);
        LOC = contains(hdr.label,'E1','IgnoreCase',1);
        ROC = contains(hdr.label,'E2','IgnoreCase',1);
        EMG1 = strcmpi(hdr.label,'Lchin');
        EMG2 = strcmpi(hdr.label,'Rchin');
        ECG2 = contains(hdr.label,'ECG2','IgnoreCase',1);
        ECG1 = contains(hdr.label,'ECG1','IgnoreCase',1);
        if sum(EMG1) > 0
            EMG = data(EMG2,:) -  data(EMG1,:);
        else
            EMG1 = contains(hdr.label,'LChinRchin','IgnoreCase',1);
            EMG = data(EMG1,:);
        end
        if sum(ECG2) > 0
            ECG = data(ECG2,:) -  data(ECG1,:);
        else
            ECG2 = contains(hdr.label,'EKG','IgnoreCase',1);
            ECG = data(ECG2,:);
        end
        EEG = data(C3,:) -  data(M2,:);
        EOGR = data(ROC,:) -  data(M1,:);
        EOGL = data(LOC,:) -  data(M2,:);
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
