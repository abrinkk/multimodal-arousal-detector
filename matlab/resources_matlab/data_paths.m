%DATA_PATHS specifies paths to data directories and lists data.
%   DATA_PATHS checks the computer profile file and selects correct data
%   directories and lists available data. 
%
%   Author: Andreas Brink-Kjaer.
%   Date: 17-Jun-2018
%

% Profile
[~,profile] = paths;

if strcmp(profile{1},'mignotserver')
    % Paths for mignot-sleep-app-01 server
    p_mros = 'E:\mros\polysomnography';
    p_mros_edf = filepath(p_mros,'edfs\visit1');
    p_mros_lab = filepath(p_mros,'annotations-events-profusion\visit1');
    p_cfs = 'E:\cfs\polysomnography';
    p_cfs_edf = filepath(p_cfs,'edfs');
    p_cfs_lab = filepath(p_cfs,'annotations-events-profusion');
    p_wsc2 = '/data2/psg/WSC_2012_2016';
    p_wsc2_edf = p_wsc2;
    p_wsc2_lab = p_wsc2;
    p_ssc = '/data2/psg/SSC/APOE';
    p_ssc_edf = p_ssc;
    p_ssc_lab = p_ssc;
    p_wsc = '/data2/psg/WSC_EDF';
    p_wsc_edf = p_wsc;
    p_wsc_lab = p_wsc;
    
    % List files
    % MrOS (empty paths)
    f_mros_edf = dir(filepath(p_mros_edf,'*.edf'));
    f_mros_edf = {f_mros_edf.name};
    f_mros_lab = dir(filepath(p_mros_lab,'*.xml'));
    f_mros_lab = {f_mros_lab.name};
    % CFS (empty paths)
    f_cfs_edf = dir(filepath(p_cfs_edf,'*.edf'));
    f_cfs_edf = {f_cfs_edf.name};
    f_cfs_lab = dir(filepath(p_cfs_lab,'*.xml'));
    f_cfs_lab = {f_cfs_lab.name};
    % WSC2
    f_wsc2_edf = dir(filepath(p_wsc2_edf,'*.edf'));
    f_wsc2_edf2 = dir(filepath(p_wsc2_edf,'*.EDF'));
    f_wsc2_edf = [{f_wsc2_edf.name} {f_wsc2_edf2.name}];
    f_wsc2_lab = dir(filepath(p_wsc2_lab,'*.csv'));
    f_wsc2_lab = {f_wsc2_lab.name};
    % Match file order
    Match_Order = cellfun(@(x) find(strcmp(cellfun(@(x) x(1:end-11),f_wsc2_lab,'Un',0),x)),cellfun(@(x) x(1:end-4),f_wsc2_edf,'Un',0),'Un',0);
    Match_Miss = cellfun(@(x) isempty(x), Match_Order);
    f_wsc2_edf(Match_Miss) = [];
    Match_Order(Match_Miss) = [];
    f_wsc2_lab = f_wsc2_lab(cell2mat(Match_Order));
    % SSC
    f_ssc_edf = dir(filepath(p_ssc_edf,'*.EDF'));
    f_ssc_edf = {f_ssc_edf.name};
    f_ssc_lab = dir(filepath(p_ssc_lab,'*.EVTS'));
    f_ssc_lab = {f_ssc_lab.name};
    % Match file order
    Match_Order = cellfun(@(x) find(strcmp(cellfun(@(x) x(1:end-5),f_ssc_lab,'Un',0),x)),cellfun(@(x) x(1:end-4),f_ssc_edf,'Un',0),'Un',0);
    Match_Miss = cellfun(@(x) isempty(x), Match_Order);
    f_ssc_edf(Match_Miss) = [];
    Match_Order(Match_Miss) = [];
    f_ssc_lab = f_ssc_lab(cell2mat(Match_Order));
    % WSC
    f_wsc_edf = dir(filepath(p_wsc_edf,'*.EDF'));
    f_wsc_edf = {f_wsc_edf.name};
    f_wsc_lab = dir(filepath(p_wsc_lab,'*.STA'));
    f_wsc_lab = {f_wsc_lab.name};
    % Match file order
    Match_Order = cellfun(@(x) find(strcmp(cellfun(@(x) x(1:end-4),f_wsc_lab,'Un',0),x)),cellfun(@(x) x(1:end-4),f_wsc_edf,'Un',0),'Un',0);
    Match_Miss = cellfun(@(x) isempty(x), Match_Order);
    f_wsc_edf(Match_Miss) = [];
    Match_Order(Match_Miss) = [];
    f_wsc_lab = f_wsc_lab(cell2mat(Match_Order));
    
    % Clear temporary variables
    clear Match_Order;
    clear Match_Miss;
else
    % Paths for local directories on external harddrive
    p_mros = 'E:\mros\polysomnography';
    p_mros_edf = filepath(p_mros,'edfs\visit1');
    p_mros_lab = filepath(p_mros,'annotations-events-profusion\visit1');
    p_cfs = 'E:\cfs\polysomnography';
    p_cfs_edf = filepath(p_cfs,'edfs');
    p_cfs_lab = filepath(p_cfs,'annotations-events-profusion');
    p_wsc2 = 'E:\wsc2\polysomnography';
    p_wsc2_edf = filepath(p_wsc2,'edfs');
    p_wsc2_lab = filepath(p_wsc2,'labels');
    p_ssc = 'E:\ssc\polysomnography';
    p_ssc_edf = filepath(p_ssc,'edfs');
    p_ssc_lab = filepath(p_ssc,'labels');
    p_wsc = 'E:\wsc\polysomnography';
    p_wsc_edf = filepath(p_wsc,'edfs');
    p_wsc_lab = filepath(p_wsc,'labels');
    
    % List files
    % MrOS
    f_mros_edf = dir(filepath(p_mros_edf,'*.edf'));
    f_mros_edf = {f_mros_edf.name};
    f_mros_lab = dir(filepath(p_mros_lab,'*.xml'));
    f_mros_lab = {f_mros_lab.name};
    % CFS
    f_cfs_edf = dir(filepath(p_cfs_edf,'*.edf'));
    f_cfs_edf = {f_cfs_edf.name};
    f_cfs_lab = dir(filepath(p_cfs_lab,'*.xml'));
    f_cfs_lab = {f_cfs_lab.name};
    % WSC2
    f_wsc2_edf = dir(filepath(p_wsc2_edf,'*.edf'));
    f_wsc2_edf = {f_wsc2_edf.name};
    f_wsc2_lab = dir(filepath(p_wsc2_lab,'*.csv'));
    f_wsc2_lab = {f_wsc2_lab.name};
    % SSC
    f_ssc_edf = dir(filepath(p_ssc_edf,'*.EDF'));
    f_ssc_edf = {f_ssc_edf.name};
    f_ssc_lab = dir(filepath(p_ssc_lab,'*.EVTS'));
    f_ssc_lab = {f_ssc_lab.name};
    % WSC
    f_wsc_edf = dir(filepath(p_wsc_edf,'*.edf'));
    f_wsc_edf = {f_wsc_edf.name};
    f_wsc_lab = dir(filepath(p_wsc_lab,'*.STA'));
    f_wsc_lab = {f_wsc_lab.name};
end
