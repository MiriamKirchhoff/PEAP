%% Load and standardize data
% Load data from original sources and save cleared data locally, only
% saving channels that are relevant for the analysis

load_data_into_standard_format


%% Split training and test data
% Split data based on time cutoff (currently: 3 min training data)
% Do this subject-wise

split_train_test_raw


%% Epoch training and test data
% Generate random epochs and epoch data accordingly
% Merge datasets of all subjects here to avoid loading time

generate_random_epochs('train')
generate_random_epochs('test')


%% Preprocess training and test data - Part 1
% Only do preprocessing steps that do not require parameter tuning, i.e. 
    % Rejection
    % Downsample
    % Detrend
    % Laplace transformation
% Save epoched data and full data array

preprocess_1('train')
preprocess_1('test')


%% Get ground truth phase estimates

get_ground_truth_phase_estimates('train', settings)
get_ground_truth_phase_estimates('test', settings)


%% Tune parameters -- Training set
% Tune parameters for all algorithms that require/allow this option
% Use first 3 min of resting state data

% 1.: PEAP
ar_opt_params = train_PEAP;

% 2.: Phastimate
% phastimate_opt_params = train_phastimate;

% 3.: ETP
ETP_opt_params = train_ETP(settings);
 

%% Filter data using all options

% Options: none, PEAP, PAR
use_filter_approaches


%% Get phase estimates using all algorithms

% Options: hilbert, phastimate, ETP, SSPE
phase_extraction_algorithms

phase_extraction_SSPE


%% Evaluate algorithms

analysis_paper
