% Script to filter data
%% Settings

path.load = '';
path.save = '';

type = 'test';


%% Load ground truth and epoched data

filepath = [path.save 'ground_truth_' type];
load([filepath '\EEG_phase_ground_truth_' type '_total'])

% load optimal parameters PEAP
filepath = [path.load '\algorithm_training_results\'];
load([filepath 'PEAP.mat'])

load('settings.mat');

data = EEG_epoch_all.data_epoched_SSPE;
time = EEG_epoch_all.t_epoched_SSPE;


%% 1. Use default filtering approach %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_current = data(end-ar_parameters_optimal.n_input+1:end, :);
time_current = EEG_epoch_all.t_epoched_SSPE(end-ar_parameters_optimal.n_input+1:end);

data_filtered_default = ...
    filtfilt(settings.filter, 1, data_current);


%% 2. Use PEAP algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Forecast AR using specified method

[data_extended_PEAP, time_extended_PEAP] = predict_ar_model( ...
    data_current, time_current, ...
    ar_parameters_optimal.order, ...
    ar_parameters_optimal.n_output , ...
    ar_parameters_optimal.approach);

%% Filter

data_filtered_PEAP = ...
    filtfilt(settings.filter, 1, data_extended_PEAP);


%% 3. Use PAR algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phast.ord = 30;
phast.hilbertwindow = 128;
phast.offset_correction = 0;
phast.edge = 65;

% Use optimized output here as well since no standard values provided in
% the paper (TODO check again)
n_interpolate = 100;
[~, ~, ~, data_all] = phastimate_adapted(...
    data_filtered_default, phast.edge, phast.ord, ...
    phast.hilbertwindow, phast.offset_correction, ...
    n_interpolate + phast.edge);

% buffer edge with nans for potential plotting
data_all_temp = [nan(size(data_all, 2), phast.edge) data_all']';

% Only use predicted data
data_ar = data_all_temp(end-n_interpolate+1:end,:);

time_extended_PAR = [time_current 0:n_interpolate-1];
% add predicted data to raw data
data_extended_PAR = [data_current; data_ar];

% timing notification
fprintf('interpolating data completed. Expired time: %.0f seconds \n', toc)


%% Filter

% get NaN trials and inf trials
nantrials = logical(sum(~isfinite(data_extended_PAR)));
data_extended_PAR(:, nantrials) = 0;
data_filtered_PAR = ...
    filtfilt(settings.filter, 1, data_extended_PAR);
data_filtered_PAR(:, nantrials) = NaN;


%% Save into array

EEG_epoch_all.data_extended_PEAP = data_extended_PEAP;
EEG_epoch_all.data_extended_PAR = data_extended_PAR;
EEG_epoch_all.data_filtered_PEAP = data_filtered_PEAP;
EEG_epoch_all.t_data_filtered_PEAP = time_extended_PEAP;
EEG_epoch_all.data_filtered_PAR = data_filtered_PAR;
EEG_epoch_all.t_data_filtered_PAR = time_extended_PAR;
EEG_epoch_all.data_filtered_default = data_filtered_default;
EEG_epoch_all.t_data_filtered_default = time_current;


%% Save common struct

% start timer
tic

% start notification
disp('Start saving common struct...')

filepath = [path.save 'filtered_' type];
% create folder if not existent
if ~exist(filepath, 'dir')
    mkdir(filepath)
end

% Save
save([filepath '\filtered_' type '_total'], ...
    "EEG_epoch_all", "n_rejection", "phase_all", "phase_full_all", '-v7.3');

% timing notification
fprintf('Saving completed. Expired time: %.0f seconds \n', toc)


