% Save prediction in raw, prob, and binary format
% Output predictions from the Multimodal Arousal Detector (MAD) in 1 second bins.
% 
% Paper: https://www.sciencedirect.com/science/article/pii/S1388245720301085.
% Code: https://github.com/abrinkk/multimodal-arousal-detector
% 
% pred_raw: Raw probabalistic predictions.
% 	line 1: arousal prediction.
% 	line 2: wake prediction.
% 	line 3: arousal predictions of data time-shifted 150 seconds.
% 	line 4: wake predictions of data time-shifted 150 seconds.
% 
% pred_prob: probabalistic predictions (center parts of 300 second segments merged from raw predictions).
% 	line 1: arousal prediction.
% 	line 2: wake prediction.
% 
% pred_binary: binary predictions (threshold probabalistic predictions).
% 	line 1: arousal predictions.
% 	line 2: wake predictions.
% 
% pred_combined: combined arousal and wake predictions (union of these).
% 	line 1: combined arousal and wake predictions.

% folders (edit for appropriate folders)
folder_raw = 'C:\Users\andre\Dropbox\MasterProject\pred\pred_raw';
folder_prob = 'C:\Users\andre\Dropbox\MasterProject\pred\pred_prob';
folder_binary = 'C:\Users\andre\Dropbox\MasterProject\pred\pred_binary';
folder_combined = 'C:\Users\andre\Dropbox\MasterProject\pred\pred_combined';

% files
f_raw = dir([folder_raw '\*.txt']);
f_raw = {f_raw.name};
n = length(f_raw);

% thresholds
T_ar = 0.225;
T_w = 0.45;

% iterate files
for i = 1:n
    % load and process data
    fname = filepath(folder_raw, f_raw{i});
    [ar_binary, w_binary] = getPred(fname,T_ar,T_w);
    [ar_prob, w_prob] = getPred(fname);
    sd = postprocess2(ar_binary,w_binary);
    
    % save data
    writematrix([ar_binary'; w_binary'],filepath(folder_binary, f_raw{i}));
    writematrix([ar_prob'; w_prob'],filepath(folder_prob, f_raw{i}));
    writematrix(sd',filepath(folder_combined, f_raw{i}));
end