function Wake = LoadWake(p_file,L,ftype)
%LOADWAKE reads sleep stages from annotation files.
%   Wake = LOADWAKE(p_file,L,ftype) loads the annotation file
%   specified in p_file and outputs wake labels.
%
%   Author: Andreas Brink-Kjaer.
%   Date: 17-Jun-2018
%
%   Input:  p_file, annotation file location
%           L, annotation length (should be based on PSG length
%           ftype, string of data source
%   Output: Wake, vector of wake sleep stages in 1 second bins.

% Determine data type
if ~exist('ftype','var')
    if contains(p_file,'mros')
        ftype = 'mros';
    elseif contains(p_file,'cfs')
        ftype = 'cfs';
    end
end

switch ftype
    case {'cfs', 'mros'}
        % Read XML file
        s = xml2struct(p_file);
        ssc = s.CMPStudyConfig.SleepStages.SleepStage;
        Wake = zeros(1,L);
        for i = 1:length(ssc)
            if contains(ssc{i}.Text,'0')
                Wake(1 + (i-1)*30:i*30) = 1;
            elseif contains(ssc{i}.Text,'1')
                Wake(1 + (i-1)*30:i*30) = 2;
            end
        end
    case 'wsc2'
        % Read csv file
        T = readtable(p_file,'FileType','text','Delimiter',';');
        T = T(:,end-1:end);
        T.Properties.VariableNames = {'Time','Event'};
        % Clear up empty event error in some Matlab versions
        T(cellfun(@isempty, T.Time),:) = [];
        % Get event times in seconds
        time_all = time2ind(T);
        Wake = zeros(1,L);
        for i = 1:size(T,1)
            if  contains(T(i,:).Event,'Stage - N2') ...
                    || contains(T(i,:).Event,'Stage - N3') || contains(T(i,:).Event,'Stage - R')
                Wake(1+round(time_all(i)):end) = 0;
            elseif contains(T(i,:).Event,'Stage - W') || contains(T(i,:).Event,'Stage - No Stage')
                Wake(1+round(time_all(i)):end) = 1;
            elseif contains(T(i,:).Event,'Stage - N1')
                Wake(1+round(time_all(i)):end) = 2;
            end
        end
    case 'ssc'
        % Read EVTS file
        [~,stageVec] = CLASS_codec.parseSSCevtsFile(p_file);
        % 30 second epoch to 1 second bins
        stageVec = repelem(stageVec,30);
        if length(stageVec) > L
            stageVec = stageVec(1:L);
        end
        Wake = ones(1,L);
        Wake(any(cell2mat(arrayfun(@(x) stageVec == x,1:5,'Un',0)),2)) = 0;
        Wake(any(stageVec == 1,2)) = 2;
    case 'wsc'
        % Read .txt file
        fid = fopen(p_file);
        STA = textscan(fid,repmat('%s',1,3));
        fclose(fid);
        stageVec = str2num(cell2mat(STA{2}));
        % 30 second epoch to 1 second bins
        stageVec = repelem(stageVec,30);
        if length(stageVec) > L
            stageVec = stageVec(1:L);
        end
        Wake = ones(1,L);
        Wake(any(cell2mat(arrayfun(@(x) stageVec == x,1:5,'Un',0)),2)) = 0;
        Wake(any(stageVec == 1,2)) = 2;
end
end