function [data, hdr] = default_edf_loader(data, hdr, channel_struct)


% Default channels
default_channel_struct = struct('EEG_Cen', {get_channel_alias([construct_channel_name({'C3'}, {'A1', 'A2', 'M1', 'M2', 'A11', 'A12', 'A21','A22','A1A2','A2A1','M1M2','M2M1', 'M11', 'M12', 'M21', 'M22', 'x', 'avg'}, {'', '-'}), {'Central'}])}, ...
    'EEG_alt', {get_channel_alias(construct_channel_name({'C4'}, {'A1', 'A2', 'M1', 'M2', 'A11', 'A12', 'A21', 'A22','A1A2','A2A1','M1M2','M2M1','M11', 'M12', 'M21', 'M22', 'x', 'avg'}, {'', '-'}))}, ...
    'REOG', {get_channel_alias(construct_channel_name({'R', 'EOG', 'E2', 'ROC', 'REOG', 'EOGR'}, {'R', 'EOG', '2', 'A1', 'A2', 'M1', 'M2', 'A11', 'A12', 'A21', 'A22', 'M11', 'M12', 'M21', 'M22', 'x', 'avg'}, {'', '-'}))}, ...
    'LEOG', {get_channel_alias(construct_channel_name({'L', 'EOG', 'E1', 'LOC', 'LEOG', 'EOGL'}, {'L', 'EOG', '2', 'A1', 'A2', 'M1', 'M2', 'A11', 'A12', 'A21', 'A22', 'M11', 'M12', 'M21', 'M22', 'x', 'avg'}, {'', '-'}))}, ...
    'EMG', {get_channel_alias(['EMG', 'Chin', construct_channel_name({'EMG', 'Chin', 'Chin1', 'Chin2', 'Chin3', 'Chin 1', 'Chin 2', 'Chin 3', 'L Chin', 'R Chin', 'LChin', 'RChin', 'Chin L', 'Chin R', 'ChinL', 'ChinR'}, {'EMG', 'Chin', 'Chin1', 'Chin2', 'Chin3', 'Chin 1', 'Chin 2', 'Chin 3', 'L Chin', 'R Chin', 'LChin', 'RChin', 'Chin L', 'Chin R', 'ChinL', 'ChinR', 'Aux12', 'Fpz'}, {'', '-'})])}, ...
    'ECG', {get_channel_alias(['ECG', 'EKG', construct_channel_name({'ECG', 'EKG', 'EKG2', 'EKG1', 'ECG2', 'ECG1', 'ECG 2', 'ECG 1', 'EKG 1', 'EKG 2', 'EKG L', 'EKG R', 'EKGR', 'EKGL', 'ECG L', 'ECG R', 'ECGR', 'ECGL', 'L EKG', 'R EKG', 'REKG', 'LEKG', 'L ECG', 'R ECG', 'RECG', 'LECG'}, {'ECG', 'EKG', 'EKG2', 'EKG1', 'ECG2', 'ECG1', 'ECG 2', 'ECG 1', 'EKG 1', 'EKG 2', 'EKG L', 'EKG R', 'EKGR', 'EKGL', 'ECG L', 'ECG R', 'ECGR', 'ECGL', 'L EKG', 'R EKG', 'REKG', 'LEKG', 'L ECG', 'R ECG', 'RECG', 'LECG'}, {'', '-'})])}, ...
    'C3', {get_channel_alias(construct_channel_name({'C3'},{'Cz','Fz','Fpz','Fp1','Fp2'},{'','-'}))}, ...
    'C4', {get_channel_alias(construct_channel_name({'C4'},{'Cz','Fz','Fpz','Fp1','Fp2'},{'','-'}))}, ...
    'ROC', {get_channel_alias(construct_channel_name({'E2 (REOG)','E2(REOG)','E2','ROC'},{'Cz','Fz','Fpz','Fp1','Fp2'},{'','-'}))}, ...
    'LOC', {get_channel_alias(construct_channel_name({'E1 (LEOG)','E1(LEOG)','E1','LOC'},{'Cz','Fz','Fpz','Fp1','Fp2'},{'','-'}))}, ...
    'A2', {get_channel_alias({'M2','A2'})}, ...
    'A1', {get_channel_alias({'M1','A1'})}, ...
    'EMG2', {get_channel_alias({'EMG 2','EMG2','EMG Aux2','EMG Chin2' 'EMG Chin 2','Chin 2','Chin2', 'Chin R', 'ChinR', 'RChin', 'R Chin'})}, ...
    'EMG1', {get_channel_alias({'EMG 1','EMG1','EMG Aux1','EMG Chin1' 'EMG Chin 1','Chin 1','Chin1', 'Chin L', 'ChinL', 'LChin', 'L Chin'})}, ...
    'ECG2', {get_channel_alias({'ECG 2', 'EKG 2', 'ECG2', 'EKG2', 'ECG I2','ECG II','ECG II2', 'EKG R', 'EKGR', 'ECG R', 'ECGR', 'R EKG', 'REKG','R ECG', 'RECG'})}, ...
    'ECG1', {get_channel_alias({'ECG 1', 'EKG 1', 'ECG1', 'EKG1','ECG I', 'EKG L', 'EKGL', 'ECG L', 'ECGL', 'L EKG', 'LEKG', 'L ECG', 'LECG'})});
        
default_field_names = fieldnames(default_channel_struct);

if ~exist('channel_struct','var')
    channel_struct = default_channel_struct;
elseif isempty(channel_struct)
    channel_struct = default_channel_struct;
elseif ~all(ismember(default_field_names, fieldnames(channel_struct)))
    for idx = 1:length(default_field_names)
        fname = default_field_names{idx};
        if ~ismember(fname, fieldnames(channel_struct))
            channel_struct.(fname) = default_channel_struct.(fname);
        end
    end
end

% EDF channel labels
file_channels = reduce_channel_name(hdr.label);

% Referneced
EEG_idx = find(ismember(file_channels,channel_struct.EEG_Cen),1,'first');
EEG_alt_idx = find(ismember(file_channels,channel_struct.EEG_alt),1,'first');
EOGR_idx = find(ismember(file_channels,channel_struct.REOG),1,'first');
EOGL_idx = find(ismember(file_channels,channel_struct.LEOG),1,'first');
EMG_idx = find(ismember(file_channels,channel_struct.EMG),1,'first');
ECG_idx = find(ismember(file_channels,channel_struct.ECG),1,'first');

% Unreferenced
C3_idx = find(ismember(file_channels,channel_struct.C3),1,'first');
C4_idx = find(ismember(file_channels,channel_struct.C4),1,'first');
ROC_idx = find(ismember(file_channels,channel_struct.ROC),1,'first');
LOC_idx = find(ismember(file_channels,channel_struct.LOC),1,'first');
A1_idx = find(ismember(file_channels,channel_struct.A1),1,'first');
A2_idx = find(ismember(file_channels,channel_struct.A2),1,'first');
EMG1_idx = find(ismember(file_channels,channel_struct.EMG1),1,'first');
EMG2_idx = find(ismember(file_channels,channel_struct.EMG2),1,'first');
ECG1_idx = find(ismember(file_channels,channel_struct.ECG1),1,'first');
ECG2_idx = find(ismember(file_channels,channel_struct.ECG2),1,'first');

% Load channels
% EEG
if ~isempty(EEG_idx)
    EEG = data{EEG_idx};
elseif ~isempty(A2_idx)
    EEG = data{C3_idx} - data{A2_idx};
    EEG_idx = C3_idx;
elseif ~isempty(C3_idx)
    EEG = data{C3_idx};
    EEG_idx = C3_idx;
elseif ~isempty(EEG_alt_idx)
    EEG = data{EEG_alt_idx};
    EEG_idx = EEG_alt_idx;
elseif ~isempty(A1_idx)
    EEG = data{C4_idx} - data{A1_idx};
    EEG_idx = C4_idx;
else
    EEG = data{C4_idx};
    EEG_idx = C4_idx;
end

% REOG
if ~isempty(EOGR_idx)
    EOGR = data{EOGR_idx};
elseif ~isempty(A2_idx)
    EOGR = data{ROC_idx} - data{A2_idx};
    EOGR_idx = ROC_idx;
else
    EOGR = data{ROC_idx};
    EOGR_idx = ROC_idx;
end

% LEOG
if ~isempty(EOGL_idx)
    EOGL = data{EOGL_idx};
elseif ~isempty(A1_idx)
    EOGL = data{LOC_idx} - data{A1_idx};
    EOGL_idx = LOC_idx;
else
    EOGL = data{LOC_idx};
    EOGL_idx = LOC_idx;
end

% Chin
if ~isempty(EMG_idx)
    EMG = data{EMG_idx};
elseif ~isempty(EMG1_idx) && ~isempty(EMG2_idx)
    EMG = data{EMG2_idx} - data{EMG1_idx};
    EMG_idx = EMG2_idx;
elseif ~isempty(EMG2_idx)
    EMG = data{EMG2_idx};
    EMG_idx = EMG2_idx;
else
    EMG = [1, -1, zeros(1, length(EEG) - 2)];
    EMG_idx = EEG_idx;
end

% ECG
if ~isempty(ECG_idx)
    ECG = data{ECG_idx};
    no_ECG = false;
elseif ~isempty(ECG2_idx) && ~isempty(ECG1_idx)
    ECG = data{ECG2_idx} - data{ECG1_idx};
    ECG_idx = ECG2_idx;
    no_ECG = false;
elseif ~isempty(ECG2_idx)
    ECG = data{ECG2_idx};
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
hdr.fs = hdr.numbersperrecord(idx) / hdr.duration; 
hdr.no_ECG = no_ECG;

end