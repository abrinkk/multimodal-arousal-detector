# multimodal-arousal-detector
This repository contains scripts for preprocessing .edf files to .txt files in the folder "matlab". The folder 'python' contains scripts for training and using the model for predicting arousals and wake.

## Getting started
Copy the repository to a local directory.

## Requirements
 * Matlab >= 2017b
 * Python >= 3.6
 * Python module numpy >= 1.14.0
 * Python module tensorflow >= 1.5.0
 
## Making new predictions
 * Locate path to folder containing .edf files of interest (p_edf) and select an output path (p_output)
 * The matlab script *LoadEDF.m* is currently able to handle .edf files from Wisconsin Sleep Cohort, Stanford Sleep Cohort, MrOS Sleep Study, and Cleveland Family Sleep Study. The script needs to be changed to handle differently structured .edf files if necessary.
 * Run "matlab" function *PreprocessNewData.m* as preprocess.PreprocessNewData(p_edf,p_output,ftype,Overwrite), where ftype can currently be one of {'mros','cfs','wsc'}, and if Overwrite = 1 files are overwritten if existing in p_output, otherwise put Overwrite = 0.
 * Run *ar_predict.py* with flags 'pathname', 'output_dir', and 'overwrite' set to input .txt file directory, output prediction directory, and binary decision to overwrite predictions, respectively.
 * In matlab, run getPred(fname,T1,T2,L) to postprocess predictions (fname: file path, T1: arousal prediction threshold = 0.225, T2: wake prediction threshold = 0.45, L: desired output length (default = prediction length)).
 
 ## Training the model
 * Edit the *profile* file to a custom name.
 * In the matlab folder, add the custom profile name in *paths.m* and specify local paths similar to other profiles.
 * In matlab/resources folder, a file called *data_paths.m* specifies directories for input data. These should be changed to list the target data.
 * The file *PreprocessData.m* in matlab/+preprocess and *LoadEDF.m* is set up to preprocess data from Wisconsin Sleep Cohort, Stanford Sleep Cohort, MrOS Sleep Study, and Cleveland Family Sleep Study. 
 * *PreprocessData.m* can now be run to preprocess data and save it in the location assigned in *paths.m*.
 * In python/ardetector *ar_config.py* the data_dir variable should be set to the same location as in *paths.m*.
 * Run *ar_train.py* with hyper-parameters specified in flags or with random hyper-parameters by setting randomgen=True in the main function. Training can be continued from pre-trained weights by setting the flag "proceed" different from 0.

