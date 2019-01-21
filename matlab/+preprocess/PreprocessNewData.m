function PreprocessNewData(p_edf,p_output,ftype,Overwrite)
%PREPROCESSNEWDATA process inputs edfs and labels
%   PREPROCESSNEWDATA inputs a set of edfs and
%   filters, resamples and normalizes data. The data is now saved in a new
%   .txt format.
%
%   Author: Andreas Brink-Kjaer.
%   Date: 15-Jan-2019
%
%   Input:  p_edf, folder locating edf files
%           p_output, output txt folder
%           ftype, file type (custom edf handling in LoadEDF.m)
%           Overwrite, to overwrite files in p_output

des_fs = 128;
dirIndex = paths;
f_edf = dir(filepath(p_edf,'*.edf'));
f_edf = {f_edf.name};
if ~exist('ftype','var')
    ftype = 'wsc';
end

parfor i = 1:length(f_edf)
    fprintf('Processsing EDFs %.0f/%.0f\n',i,length(f_edf));
    try
        if ~Overwrite && exist(filepath(p_output,[f_edf{i}(1:end-4) '.txt']),'file')
            continue;
        end
        % Load and preprocess data
        [hdr,data] = LoadEDF(filepath(p_edf,f_edf{i}),ftype);
        [hdr,data] = preprocess.resampleData(data,hdr,des_fs);
        data = preprocess.filterData(data,hdr,1);
        
        % Load labels
        ar_seq = zeros(1,size(data,2)/des_fs);
        W = ar_seq;
        
        % Save
        exportData(data,ar_seq,W,hdr,filepath(p_output,[f_edf{i}(1:end-4) '.txt']));
    catch me
        disp(me.message);
    end
end
end
