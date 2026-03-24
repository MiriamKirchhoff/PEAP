function generate_random_epochs(type)


if nargin < 1
    type = 'train';
end

%% Set paths

path.load = '';
path.save = '';


%% Set epoch windows

win.epoch       =       [-0.705 0.00001];
win.epoch_SSPE  =       [-2.065 0.00001];
win.epoch_validation =  [-.8 .3];


%% Set loading options

subjects = 1:50;

datasets = {'young_control', 'stroke_ipsi', 'stroke_contra', 'elderly_control'};


%% Initialization

eeglab

% iterate over datasets
for idx_dataset = 1:length(datasets)

    dataset = datasets{idx_dataset};
    fprintf('\n \n START DATASET: %s \n \n', dataset);

    % iterate over subjects
    for idx_subject = subjects

        formatSpec = '%03.0f';


        %% Load subset resting state data

        % start timer
        tic

        % start notification
        fprintf('\nStarting participant %03.0f. \n', idx_subject)

        filepath = [path.load 'raw_' type '\' dataset '\EEG_rs_' type '_' dataset ...
            '_'  num2str(idx_subject,formatSpec) '.mat'];

        if ~exist(filepath, 'file')
            fprintf('\nNo data found for participant %03.0f. \n', idx_subject)
            continue
        end

        % load EEG data
        disp('Start loading...')
        load(filepath);

        % timing notification
        fprintf('loading completed. Expired time: %.0f seconds \n', toc)


        %% Create events for epoching

        isi.mean = 1.5;
        isi.jitter = 0.25;
        % maximum amount of epochs possible with dataset size
        isi.n_possible = round(EEG.xmax/(isi.mean - isi.jitter));

        % get random inter-stimulus intervals between 0 and 1
        isi.all = rand(1, isi.n_possible);
        % get mean 0 and range [-1 1], then scale to +- jitter
        isi.all = (isi.all*2 - 1) * isi.jitter;
        % add mean isi
        isi.all = isi.all + isi.mean;

        isi.eventtimes = round(cumsum(isi.all)*EEG.srate);

        % drop all events that are beyond the measurement
        isi.min_timepoint = -win.epoch_SSPE(1)*EEG.srate;
        isi.max_timepoint = EEG.pnts - win.epoch_validation(end)*EEG.srate;

        isi.exclude = isi.eventtimes < isi.min_timepoint | ...
            isi.eventtimes > isi.max_timepoint;
        isi.all(isi.exclude) = [];
        isi.eventtimes(isi.exclude) = [];


        %% Include events in EEGlab struct

        event_type = ['random, ISI ' num2str(isi.mean) ' , jitter ' num2str(isi.jitter)];
        n_events = numel(isi.eventtimes);
        EEG.event = struct( ...
            'type', repmat({event_type}, 1, n_events), ...
            'latency', num2cell(isi.eventtimes), ...
            'urevent', num2cell(1:n_events) ...
            );

        EEG.urevent = struct( ...
            'type', repmat({event_type}, 1, n_events), ...
            'latency', num2cell(isi.eventtimes) ...
            );

        EEG = eeg_checkset(EEG, 'eventconsistency');
        EEG = eeg_checkset(EEG, 'makeur');


        %% Epoch based on random ISIs

        % cut off beyond 0 ms first since EEGlab can not cut off before event
        EEG_epoched = pop_epoch(EEG, event_type, win.epoch);
        %EEG_epoched = pop_select(EEG_epoched, 'time', win.epoch);

        % Generate long epochs
        EEG_epoched_SSPE = pop_epoch(EEG, event_type, win.epoch_SSPE);
        %EEG_epoched_SSPE = pop_select(EEG_epoched_SSPE, 'time', win.epoch_SSPE);

        % Generate ground truth raw data for later validation
        EEG_epoched_raw = pop_epoch(EEG, event_type, win.epoch_validation);


        %% Reject noisy trials

        outlier_threshold = 150;

        % start timer
        tic

        % start notification
        disp('Start rejecting trials...')

        % calculate range of all channels across time
        ranges = squeeze(range(EEG_epoched.data,2));

        % consider only highest idx
        ranges = max(ranges);

        % outlier detection
        outlier_location = find(ranges > outlier_threshold);
        n_rejection = length(outlier_location);

        if n_rejection > 0
            % reject epochs
            keep_idx = find(ranges <= outlier_threshold);
            EEG_epoched = pop_select(EEG_epoched, 'trial', keep_idx);
            EEG.event(outlier_location) = [];
            EEG_epoched_SSPE = pop_select(EEG_epoched_SSPE, 'trial', keep_idx);
            EEG_epoched_raw = pop_select(EEG_epoched_raw, 'trial', keep_idx);
        end

        % timing notification
        fprintf('Rejection completed. Trials: %.0f. Expired time: %.0f seconds \n', n_rejection, toc)


        %% Save data

        filepath = [path.save 'raw_epoched_' type '\' dataset];
        % create folder if not existent
        if ~exist(filepath, 'dir')
            mkdir(filepath)
        end

        % Save
        save([filepath '\EEG_epoched' '_' type '_' dataset '_'  num2str(idx_subject,formatSpec)], ...
            "EEG", "EEG_epoched", "EEG_epoched_SSPE", "EEG_epoched_raw", "n_rejection", '-v7.3');

    end % idx_subject = subjects  % iterate over subjects
end % for idx_dataset = 1:length(datasets)

disp('Finished generating random epochs')