# multimodal-arousal-detector
This repository contains scripts for preprocessing .edf files to .txt files in the folder 'matlab'. The folder 'python' and 'python_tf2' contains scripts for training and using the model for predicting arousals and wake.

Please see the [publication](https://www.sciencedirect.com/science/article/pii/S1388245720301085).

## Getting started
 * Copy the repository to a local directory.
 * Create a folder to store preprocessed .edf files.
 * Create a folder to store output arousal and wake predictions.
 
## Requirements (Updated code to work with Tensorflow v. 2)
 * Matlab >= 2017b
 * Python = 3.8 
 * Python module Tensorflow v. 2.10.0

## Making new predictions (Updated code to work with Tensorflow v. 2)
 * Locate path to folder containing .edf files of interest (p_edf) and select an output path (p_output)
 * The matlab script *LoadEDF.m* should be able to recognize all channel names in English (e.g. 'C3-A2' or 'LEOG') and can handle missing ECG data. If not you can make a new 'ftype' case to handle recognizing special channel names.
 * Run "matlab" function *PreprocessNewData.m* as `preprocess.PreprocessNewData(p_edf, p_output, ftype, Overwrite)`, where 'ftype' should be set to 'default' or something else for a specific case, and if 'Overwrite' = 1 files are overwritten if existing in p_output, otherwise put 'Overwrite' = 0.
 * Run *ar_predict.py* with flags 'pathname', 'output_dir', and 'overwrite' set to input .txt file directory, output prediction directory, and the binary decision to overwrite predictions, respectively. 'pathname' should be the same as 'p_output' when *PreprocessNewData.m* was run.
 * In matlab, run `getPred(fname,T1,T2,L)` to postprocess predictions (fname: file path, T1: arousal prediction threshold = 0.225, T2: wake prediction threshold = 0.45, L: desired output length (default = prediction length)).
 * The matlab script `save_prediction_format` can save the predictions in various formats, which can be easier to work with subsequently.

## Requirements (Original version with Tensorflow v. 1.5.0)
 * Matlab >= 2017b
 * Python >= 3.6 (see requirements.txt for conda environemnt list)
 * Python module numpy >= 1.14.0
 * Python module tensorflow v. 1.5.0 (currently incompatible with version 2.0)
   - (Tensorflow v. 1.5.0 requires cuDNN 7.0.4 and CUDA 9.0)

## Making new predictions (Original version with Tensorflow v. 1.5.0)
 * Locate path to folder containing .edf files of interest (p_edf) and select an output path (p_output)
 * The matlab script *LoadEDF.m* is currently able to handle .edf files from Wisconsin Sleep Cohort, Stanford Sleep Cohort, MrOS Sleep Study, and Cleveland Family Sleep Study. The script needs to be changed to handle differently structured .edf files if necessary.
 * Run "matlab" function *PreprocessNewData.m* as `preprocess.PreprocessNewData(p_edf,p_output,ftype,Overwrite)`, where ftype can currently be one of {'mros','cfs','wsc'}, and if Overwrite = 1 files are overwritten if existing in p_output, otherwise put Overwrite = 0.
 * Run *ar_predict.py* with flags 'pathname', 'output_dir', and 'overwrite' set to input .txt file directory, output prediction directory, and binary decision to overwrite predictions, respectively. 'pathname' should be the same as 'p_output' when *PreprocessNewData.m* was run.
 * In matlab, run `getPred(fname,T1,T2,L)` to postprocess predictions (fname: file path, T1: arousal prediction threshold = 0.225, T2: wake prediction threshold = 0.45, L: desired output length (default = prediction length)).
 * The matlab script `save_prediction_format` can save the predictions in various formats, which can be easier to work with subsequently.
 
## Data access
Polysomnography data from the MrOS Sleep Study and Cleveland Family Sleep Study is available upon request from the National Sleep Research Resource (NSRR).

 ## Training the model
 * In matlab/resources folder, a file called *data_paths.m* specifies directories for input data. These should be changed to list the target data.
 * The file *PreprocessData.m* in matlab/+preprocess and *LoadEDF.m* is set up to preprocess data from Wisconsin Sleep Cohort, Stanford Sleep Cohort, MrOS Sleep Study, and Cleveland Family Sleep Study. 
 * *PreprocessData.m* can now be run to preprocess data and save it in the location assigned in *paths.m*.
 * In python/ardetector *ar_config.py* the data_dir variable should be set to the same location as in *paths.m*.
 * Run *ar_train.py* with hyper-parameters specified in flags or with random hyper-parameters by setting randomgen=True in the main function. Training can be continued from pre-trained weights by setting the flag "proceed" different from 0.
