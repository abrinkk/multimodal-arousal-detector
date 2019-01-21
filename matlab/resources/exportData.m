function exportData(data,ar,W,hdr,p_file)
%EXPORTDATA writes a .txt file with processed data.
%   EXPORTDATA(data,ar,W,hdr,p_file) inputs resampled and filtered EEG, EOG
%   and EMG data. The data is reshaped and saved in a .txt format suited to
%   pass to the ar_reader.py Python script used to feed data to the model.
%
%   Author: Andreas Brink-Kjaer.
%   Date: 17-Jun-2018
%
%   Input:  data, resampled and filtered EGE, EOG and EMG data
%           ar, arousal annotations (if not available, then pass zeros)
%           W, wake annotatinos (if not available, then pass zeros)
%           hdr, .edf header
%           p_file, destination file

% File identifier
FID = fopen(p_file,'w');
% Data reshape
data = [reshape(data(1,:),hdr.fs(1),[])' reshape(data(2,:),hdr.fs(2),[])' ...
    reshape(data(3,:),hdr.fs(3),[])' reshape(data(4,:),hdr.fs(4),[])' ar' W'];
% Write data to file
fprintf(FID,[repmat('%.2f,',1,size(data,2)-2) '%.0f, %.0f\n'],data');
% Close file
fclose(FID);
end