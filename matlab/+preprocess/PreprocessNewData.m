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
if strcmp(ftype,'dreem')
    f_edf = dir(filepath(p_edf,'*.h5'));
    f_edf = {f_edf.name};
    f_edf = unique(f_edf);
else
    f_edf_1 = dir(filepath(p_edf,'*.edf'));
    f_edf_2 = dir(filepath(p_edf,'*.EDF'));
    f_edf = [f_edf_1; f_edf_2];
    f_edf = {f_edf.name};
    f_edf = unique(f_edf);
end
if ~exist('ftype','var')
    ftype = 'wsc';
end

% parpool(2);
for i = 1:length(f_edf)
    fprintf('Processsing EDFs %.0f/%.0f\n',i,length(f_edf));
    try
        f_edf_i = f_edf{i};
        if strcmp(ftype, 'dreem')
            f_edf_i_short = f_edf_i(1:end-3);
        else
            f_edf_i_short = f_edf_i(1:end-4);
        end
        if ~Overwrite && exist(filepath(p_output,[f_edf_i_short '.txt']),'file')
            continue;
        end
        % Load and preprocess data
        [hdr,data] = LoadEDF(filepath(p_edf,f_edf_i),ftype);
        data = data(:,1:max(hdr.fs)*floor(size(data,2)/max(hdr.fs)));
        % Does the data contain ECG
        if isfield(hdr, 'no_ECG')
            filter_EKG = ~hdr.no_ECG;
        else
            filter_EKG = 1;
        end
        [hdr,data] = preprocess.resampleData(data,hdr,des_fs);
        data = preprocess.filterData(data,hdr,filter_EKG);
        % Load labels
        ar_seq = zeros(1,size(data,2)/des_fs);
        W = ar_seq;
        
        % Save
        exportData(data,ar_seq,W,hdr,filepath(p_output,[f_edf_i_short '.txt']));
    catch me
        disp(me.message);
    end
end
end
