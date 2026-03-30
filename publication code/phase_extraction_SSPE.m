% This code is retrieved and altered from the Wodeyar (2021) publication 

%% Clear

clear all

%% Settings - general

path.load = '';
path.save = path.load;

type = 'test';
target_time = -1;     % [ms]


%% Load results of phase extraction of other algorithms

filepath = [path.load 'phase_' type];
load([filepath '\phase_' type '_total'])


%% SSPE settings

% setting up initial parameters to start the causalPhaseEM code
% setting up to estimate three oscillators - helps to account for the harmonic of the 11 Hz
SSPE_settings.freqs = [2,11,22]; 
SSPE_settings.Fs = 1000;

% in a pinch this can be initialized to 0.99 to start
SSPE_settings.ampVec = [.99,.99,.99]; 

% its important to use a value of this that is in the same ballpark scale
SSPE_settings.sigmaFreqs = [.01,.01,.01]; 
SSPE_settings.sigmaObs = 1;
SSPE_settings.window = 1000;

% adapted from [8 14] to match other analyses
SSPE_settings.lowFreqBand = [9, 13]; 

%% Calculate phase based on SSPE

% Based on unfiltered data

phase_full_all.t_SSPE = EEG_epoch_all.t_epoched_SSPE;

prev_msg = '';
n_trials = size(EEG_epoch_all.data_epoched_SSPE,2);
for idx_trial = 1:n_trials
    % For overview only
    msg = sprintf('Processing trial %d of %d \n', idx_trial, n_trials);
    fprintf([repmat('\b', 1, length(prev_msg)), msg]);
    prev_msg = msg;

    chunk = EEG_epoch_all.data_epoched_SSPE(:,idx_trial);
    [phase_chunk,phaseBounds,allX_full,phaseWidth,returnParams] = causalPhaseEM_MKmdl_noSeg(chunk, SSPE_settings,0);
    phase_all.SSPE(idx_trial) = phase_chunk(end);
    phase_full_all.SSPE(idx_trial,:) = phase_chunk;
end

fprintf('\n');
clear chunk msg prev_msg

% plot_KF(y,initParams)


%% Save common struct

% start timer
tic

% start notification
disp('Start saving common struct...')

filepath = [path.save 'phase_' type];
% create folder if not existent
if ~exist(filepath, 'dir')
    mkdir(filepath)
end

% Save
save([filepath '\phase_' type '_total'], ...
    "EEG_epoch_all", "n_rejection", "phase_all", "phase_full_all", '-v7.3');

% timing notification
fprintf('Saving completed. Expired time: %.0f seconds \n', toc)
