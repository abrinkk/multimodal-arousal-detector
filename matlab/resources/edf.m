classdef edf < handle
    %EDF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dataSource
        dataRaw
        dataFiltered
        hdr
        loaded
    end
    
    properties
        Channel_EEG_Cen = {'C3M2' 'C3M1' 'C4M2' 'C4M1' 'C4AVG' 'C3AVG','C3A2','C4A1','C3A1','C4A2','C3x','C4x'}
        Channel_EEG_Fro = {'F3M2' 'F3M1' 'F3AVG' 'F4AVG' 'F4M1','F1A2','F2C4','FzA2','FzA1','F3x','F4x'}
        Channel_EEG_Occ = {'O1M2' 'O1M1' 'O2M2' 'O2M1' 'O2M1','O1AVG','O1x','O2x'}
        Channel_EMG = {'Chin1Chin2' 'Chin1Chin3' 'Chin3Chin2' 'ChinEMG'}
        Channel_EKG = {'ECG' 'EKG1AVG' 'EKG1EKG2' 'LLEG1EKG2' 'EKG'}
        Channel_LEOG = {'LEOGM2','LOCM1','LOCA2','LOCA1','LEOGx'}
        Channel_REOG = {'REOGM1','ROCM1','ROCA1','ROCA2','REOGx'}
    end
    
    methods
        
        function [ obj ] = edf( source )
            obj.dataSource = source;
            obj.dataRaw = [];
            obj.dataFiltered = [];
            obj.hdr = [];
            obj.loaded = false;
            obj.loadHeader(edfread(obj.dataSource));
        end
        
        % Function to load everything and resample for APOE data
        function [ x ] = resample_get_all(obj, fro, raw )
            if nargin < 2
                fro = 1;
            end
            if nargin < 3
                raw = false;
            end
            if fro == 1
                channels = {obj.Channel_EEG_Cen, obj.Channel_EEG_Fro, ...
                    obj.Channel_EEG_Occ, obj.Channel_EMG, ...
                    obj.Channel_LEOG, obj.Channel_REOG, obj.Channel_EKG};
            else
                channels = {obj.Channel_EEG_Cen, ...
                    obj.Channel_EEG_Occ, obj.Channel_EMG, ...
                    obj.Channel_LEOG, obj.Channel_REOG, obj.Channel_EKG};
            end
            obj.loadData();
            % Fill in zeros to match desired fs
            des_fs = 200;
            max_fs = max(obj.hdr.fs);
            if des_fs > max_fs
                obj.dataRaw = [obj.dataRaw zeros(size(obj.dataRaw,1),size(obj.dataRaw,2)*(des_fs/max_fs - 1))];
                obj.dataFiltered = nan(size(obj.dataRaw));
                max_fs = des_fs;
            end
            x = zeros(length(channels), size(obj.dataRaw,2));
            fs = zeros(length(channels),1);
            idxc = fs;
            for i = 1:length(channels)
                idx = ismember(obj.hdr.label, channels{i});
                assert(any(idx));
                if sum(idx) > 1
                    % if there is more than one match (e.g. nasal pressure
                    % AND airflow can exist at same time, even though in
                    % same category), then select the one first listed
                    keys = find(ismember(channels{i}, obj.hdr.label));
                    key = channels{i}(keys(1));
                    idx = ismember(obj.hdr.label, key);
                end
                idxc(i) = find(idx);
                fs(i) = obj.hdr.fs(idx);
                x(i,:) = obj.dataRaw(idx,:);
            end
            % Resampling assuming zero padding at end matching
            % highest sampling frequency
            obj.hdr.fs = des_fs*ones(size(obj.hdr.fs));
            for i = 1:length(channels)
                dataEnd = round(fs(i)/max_fs*length(x));
                dataClip = x(i,1:dataEnd);
                % Check up/down sample factor
                [p,q] = rat(des_fs/fs(i));
                if fs(i) ~= des_fs
                    dataResample = resample(dataClip,p,q);
                    x(i,1:round(dataEnd*p/q)) = dataResample;
                else
                    x(i,1:dataEnd) = dataClip;
                end
                obj.dataRaw(idxc(i),:) = x(i,:);
            end
            % Filtering
            for i = 1:length(channels)-1
                if ~raw
                    if all(isnan(obj.dataFiltered(idxc(i),:)))
                        obj.dataFiltered(idxc(i),:) = obj.filterRow(channels{i}, x(i,:));
                    end
                    x(i,:) = obj.dataFiltered(idxc(i),:);
                end
            end
            dataEnd = round(des_fs/max_fs*length(x));
            x = x(:,1:dataEnd);
            obj.dataRaw = obj.dataRaw(:,1:dataEnd);
            obj.dataFiltered = obj.dataFiltered(:,1:dataEnd);
        end
        
        function [ x ] = get( obj, channels, raw )
            if nargin < 3
                raw = false;
            end
            obj.loadData();
            x = zeros(length(channels), size(obj.dataRaw,2));
            for i = 1:length(channels)
                idx = ismember(obj.hdr.label, channels{i});
                assert(any(idx));
                if sum(idx) > 1
                    % if there is more than one match (e.g. nasal pressure
                    % AND airflow can exist at same time, even though in
                    % same category), then select the one first listed
                    keys = find(ismember(channels{i}, obj.hdr.label));
                    key = channels{i}(keys(1));
                    idx = ismember(obj.hdr.label, key);
                end
                x(i,:) = obj.dataRaw(idx,:);
                if ~raw
                    if all(isnan(obj.dataFiltered(idx,:)))
                        obj.dataFiltered(idx,:) = obj.filterRow(channels{i}, x(i,:));
                    end
                    x(i,:) = obj.dataFiltered(idx,:);
                end
            end
        end
        
        function [ sampleRate ] = fs( obj, index )
            sampleRate = mode(obj.hdr.fs);
        end
        
    end
    
    methods
        
        function [ x ] = filterRow( obj, id, x )
            % EEG
            if any(ismember(obj.Channel_EEG_Cen, id)) || ...
                    any(ismember(obj.Channel_EEG_Fro, id)) || ...
                    any(ismember(obj.Channel_LEOG, id)) || ...
                    any(ismember(obj.Channel_REOG, id))
                x = obj.filterEKGArtefact(x);
                x = obj.filterEEGbp(x);
                
                % EMG Chin
            elseif any(ismember(obj.Channel_EMG, id))
                x = obj.filterEKGArtefact(x);
                x = obj.filterEMGbp(x);
                
            end
        end
        
        function y = filterEMGbp( obj, x )
            b = obj.emgFilter(obj.fs);
            y = filtfilt(b, 1, x);
        end
        
        function y = filterEEGbp( obj, x )
            b = obj.eegFilter(obj.fs);
            y = filtfilt(b, 1, x);
        end
        
        function y = filterEKGArtefact( obj, x )
            sig_data = x;
            ref_data = obj.get({obj.Channel_EKG}, 1);
            least = min([length(sig_data) length(ref_data)]);
            e = filter.anc_rls_m(sig_data(1:least)', ref_data(1:least)');
            % Can contain nan for flat segments
            e(isnan(e)) = sig_data(isnan(e));
            y = zeros(size(sig_data));
            y(1:length(e)) = e;
            y = y(:)';
        end
        
        function y = filterElectrodePopping( obj, x )
            y = wsc.filter.electrode(x);
        end
        
    end
    
    methods (Access=protected)
        
        function loadData( obj )
            if ~obj.loaded
                [hdr_, data_] = edfread(obj.dataSource);
                obj.dataRaw = data_;
                obj.dataFiltered = nan(size(data_));
                obj.loadHeader(hdr_);
                obj.loaded = 1;
            end
        end
        
        function loadHeader( obj, hdr_ )
            obj.hdr = hdr_;
            if ~ismember('fs',fieldnames(obj.hdr))
                obj.hdr.fs = obj.hdr.samples ./ obj.hdr.duration;
            end
        end
        
    end
    
    methods (Access=public)
        
        function assert( obj )
            % Assert all required channels exist
            obj.loadData;
            fnames = fieldnames(obj);
            cnames = find(~cellfun(@isempty,regexpi(fnames,'Channel_(.*)')));
            for i = 1:length(cnames)
                key = obj.(fnames{cnames(i)});
                message = sprintf('CHANNEL_MISSING: %s', fnames{cnames(i)});
                assert(any(ismember(key, obj.hdr.label)), message);
            end
            
            % Assert there is not too many artifacts
            artifacts = obj.getArtifacts();
            artifactRatio = mean(any(artifacts,1));
            message = sprintf('NOISE_TOO_HIGH: %.2f', artifactRatio);
            assert(artifactRatio < 0.75, message);
        end
        
        function [ h ] = genOverviewPlot( obj )
            h = figure;
            h.Position(3:4) = [700 900];
            centerfig(h);
            
            fnames = fieldnames(obj);
            cnames = find(~cellfun(@isempty,regexpi(fnames,'Channel_(.*)')));
            for i = 1:length(cnames)
                key = obj.(fnames{cnames(i)});
                [x] = obj.get({key});
                [fx,fy] = fftcustom(x, obj.fs);
                subplot(length(cnames)+1, 2, 1+(i-1)*2);
                plot(x);
                box on, grid minor, axis tight;
                ylabel(fnames{cnames(i)},'Interpreter','None');
                subplot(length(cnames)+1, 2, i*2);
                plot(fx,fy);
                box on, grid minor, axis tight;
            end
            
            % Artifacts and arousals
            subplot(length(cnames)+1, 2, length(cnames)*2 + [1 2]);
            arousals = obj.getArousals();
            artifacts = obj.getArtifacts();
            plot(1:length(arousals),arousals,1:length(artifacts),any(artifacts,1));
            box on, grid minor, axis tight;
            legend('Arousals','Artifacts');
            ylabel('Misc.');
        end
        
        function [ n ] = channelNoise( obj, key, x )
            n = false(size(x));
            
            if ismember(key, {'Channel_EEG' 'Channel_EMG' 'Channel_LEG'})
                srcSig = x';
                params = struct;
                params.samplerate = obj.fs;
                params.win_length_sec = 3;
                params.win_interval_sec = 3;
                stageStruct = [];
                eventStruct = detection.detection_artifact_electrode_pop(srcSig, params, stageStruct);
                indexed = obj.indexEventStruct(eventStruct);
                n = any([n; indexed], 1);
                
            end
        end
        
        function [ noise ] = getArtifacts( obj )
            fnames = fieldnames(obj);
            cnames = find(~cellfun(@isempty,regexpi(fnames,'Channel_(.*)')));
            noise = false(size(obj.dataRaw));
            for i = 1:length(cnames)
                key = obj.(fnames{cnames(i)});
                x = obj.get({ key });
                n = obj.channelNoise(fnames{cnames(i)}, x);
                noise(i,:) = n;
            end
        end
        
    end
    
    methods (Static)
        
        function [ indexed ] = indexEventStruct( eventStruct )
            indexed = false(1, size(eventStruct.new_data,1));
            for i = 1:size(eventStruct.new_events,1)
                i1 = eventStruct.new_events(i,1);
                i2 = eventStruct.new_events(i,2);
                if ~any(isnan([i1 i2])) && ~any([i1 i2] > length(indexed))
                    indexed(i1:i2) = true;
                end
            end
        end
        
        function b = emgFilter(fs)
            Fstop1 = 5;               % First Stopband Frequency
            Fpass1 = 10;              % First Passband Frequency
            Fpass2 = 45;              % Second Passband Frequency
            Fstop2 = 50;              % Second Stopband Frequency
            Dstop1 = 0.01;            % First Stopband Attenuation
            Dpass  = 0.057501127785;  % Passband Ripple
            Dstop2 = 1e-05;           % Second Stopband Attenuation
            dens   = 20;              % Density Factor
            [N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(fs/2), [0 1 0], [Dstop1 Dpass Dstop2]);
            b  = firpm(N, Fo, Ao, W, {dens});
        end
        
        function b = eegFilter(fs)
            Fpass = 45;              % Passband Frequency
            Fstop = 50;              % Stopband Frequency
            Dpass = 0.057501127785;  % Passband Ripple
            Dstop = 1e-05;           % Stopband Attenuation
            dens  = 20;              % Density Factor
            [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(fs/2), [1 0], [Dpass, Dstop]);
            b  = firpm(N, Fo, Ao, W, {dens});
        end
        
    end
    
end

