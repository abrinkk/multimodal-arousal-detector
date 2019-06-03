function SSC = LoadSSC(p_file,L,ftype)
%LOADSSC reads sleep stages from annotation files.
%   SSC = LOADSSC(p_file,L,ftype) loads the annotation file
%   specified in p_file and outputs the sleep stages.
%
%   Author: Andreas Brink-Kjaer.
%   Date: 17-Jun-2018
%
%   Input:  p_file, annotation file location
%           L, annotation length (should be based on PSG length
%           ftype, string of data source
%   Output: SSC, vector of sleep stages in 1 second bins.

% Unknown length
L_unk = 0;
if isempty(L)
    L = 20*60*60;
    L_unk = 1;
end

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
        SSC = cellfun(@(x) str2num(x.Text), ssc);
        % 30 second epochs to 1 second bins
        SSC = repelem(SSC,30);
    case 'wsc2'
        % Read csv file
        T = readtable(p_file,'FileType','text','Delimiter',';');
        T = T(:,end-1:end);
        T.Properties.VariableNames = {'Time','Event'};
        % Clear up empty event error in some Matlab versions
        T(cellfun(@isempty, T.Time),:) = [];
        % Get event times in seconds
        time_all = time2ind(T);
        if L_unk == 1
            L = ceil(time_all(end));
        end
        SSC = zeros(1,L);
        for i = 1:size(T,1)
            if contains(T(i,:).Event,'Stage - W') || contains(T(i,:).Event,'Stage - No Stage')
                SSC(1+round(time_all(i)):end) = 0;
            elseif contains(T(i,:).Event,'Stage - N1')
                SSC(1+round(time_all(i)):end) = 1;
            elseif contains(T(i,:).Event,'Stage - N2')
                SSC(1+round(time_all(i)):end) = 2;
            elseif contains(T(i,:).Event,'Stage - N3')
                SSC(1+round(time_all(i)):end) = 3;
            elseif contains(T(i,:).Event,'Stage - R')
                SSC(1+round(time_all(i)):end) = 5;
            end
        end
    case 'ssc'
        % Read .EVTS file
        [~,stageVec] = CLASS_codec.parseSSCevtsFile(p_file);
        % 30 second epochs to 1 second bins
        stageVec = repelem(stageVec,30);
        if length(stageVec) > L
            stageVec = stageVec(1:L);
        end
        SSC = stageVec;
        SSC(any(stageVec == 6:7,2)) = 0;
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
        SSC = stageVec;
        SSC(~any(cell2mat(arrayfun(@(x) stageVec == x,1:5,'Un',0)),2)) = 0;
end
end
