function preprocess_1(type)


if nargin < 1
    type = 'train';
end

%% Set paths

path.load = '';
path.save = '';


%% Set loading options

subjects = 1:50;

datasets = {'young_control', 'stroke_ipsi', 'stroke_contra', 'elderly_control'};


%% Initialization

eeglab


%% Iterate over subjects


for idx_dataset = 1:length(datasets)

    dataset = datasets{idx_dataset};
    fprintf('\n \n START DATASET: %s \n \n', dataset);

    for idx_subject = subjects  % iterate over subjects

        formatSpec = '%03.0f';


        %% Load standardized resting state REFTEP data

        % start timer
        tic

        % start notification
        fprintf('\nStarting participant %03.0f. \n', idx_subject)


        filepath = [path.load 'raw_epoched_' type '\' dataset '\EEG_epoched' '_' type '_' dataset '_'  ...
            num2str(idx_subject,formatSpec) '.mat'];
        if ~exist(filepath, 'file')
            fprintf('\nNo data found for participant %03.0f. \n', idx_subject)
            continue
        end

        n_rejection = [];
        load(filepath);

        % timing notification
        fprintf('loading completed. Expired time: %.0f seconds \n', toc)

        %% Settings for montage

        switch dataset
            case 'stroke_ipsi'
                % Ipsilaterial hemisphere
                ipsi_loc = {'C4' 'C4' 'C4' 'C4' 'C3' 'C3' 'C3' 'C3' '' 'C4'};
                hjorth.channel = ipsi_loc{idx_subject};
                clear ipsi_loc
            case 'stroke_contra'
                % Contralateral hemisphere
                contra_loc = {'C3' 'C3' 'C3' 'C3' 'C4' 'C4' 'C4' 'C4' '' 'C3'};
                hjorth.channel = contra_loc{idx_subject};
                clear contra_loc
            otherwise
                hjorth.channel = 'C3';
        end

        switch hjorth.channel  % get surrounding variables for hjorth filter
            case 'C3'
                hjorth.electrodes = {'FCC5h', 'FCC3h', 'CCP5h', 'CCP3h'};
            case 'C4'
                hjorth.electrodes = {'FCC6h', 'FCC4h', 'CCP6h', 'CCP4h'};
        end


        %% Do initial preprocessing: Downsample, detrend, laplacian montage

        EEG                 = prep_ds_dt_lap(EEG, hjorth);
        EEG_epoched         = prep_ds_dt_lap(EEG_epoched, hjorth);
        EEG_epoched_raw     = prep_ds_dt_lap(EEG_epoched_raw, hjorth);
        EEG_epoched_SSPE    = prep_ds_dt_lap(EEG_epoched_SSPE, hjorth);


        %% Save final data

        % start timer
        tic

        % start notification
        disp('Start saving data...')

        filepath = [path.save 'prep_1_' type '\' dataset];
        % create folder if not existent
        if ~exist(filepath, 'dir')
            mkdir(filepath)
        end

        % Save
        save([filepath '\EEG_prep_1_' type '_' dataset '_'  num2str(idx_subject,formatSpec)], ...
            "EEG", "EEG_epoched", "EEG_epoched_SSPE", "EEG_epoched_raw", "n_rejection", '-v7.3');

        % timing notification
        fprintf('Saving completed. Expired time: %.0f seconds \n', toc)


        %% Merge into common struct 

        % initialize
        if ~exist('EEG_epoch_all', 'var')
            EEG_epoch_all.t_epoched     = EEG_epoched.times;
            EEG_epoch_all.t_epoched_raw = EEG_epoched_raw.times;
            EEG_epoch_all.t_epoched_SSPE= EEG_epoched_SSPE.times;
            EEG_epoch_all.data_epoched     = [];
            EEG_epoch_all.data_epoched_raw = [];
            EEG_epoch_all.data_epoched_SSPE= [];
            EEG_epoch_all.t_cont            = {};
            EEG_epoch_all.data_cont         = {};
            EEG_epoch_all.events_cont        = {};
            EEG_epoch_all.subject           = [];
            EEG_epoch_all.dataset           = [];
            EEG_epoch_all.fs                = EEG.srate;
            EEG_epoch_all.n_rejection       = [];
        end

        EEG_epoch_all.data_epoched     = [EEG_epoch_all.data_epoched, squeeze(EEG_epoched.data)];
        EEG_epoch_all.data_epoched_raw = [EEG_epoch_all.data_epoched_raw, squeeze(EEG_epoched_raw.data)];
        EEG_epoch_all.data_epoched_SSPE= [EEG_epoch_all.data_epoched_SSPE, squeeze(EEG_epoched_SSPE.data)];
        EEG_epoch_all.t_cont           = [EEG_epoch_all.t_cont; {EEG.times}];
        EEG_epoch_all.data_cont        = [EEG_epoch_all.data_cont; {EEG.data}];
        EEG_epoch_all.events_cont      = [EEG_epoch_all.events_cont; [EEG.event.latency]];
        EEG_epoch_all.subject          = [EEG_epoch_all.subject, ones(1,EEG_epoched.trials)*idx_subject + 100*idx_dataset];
        EEG_epoch_all.n_rejection      = [EEG_epoch_all.n_rejection; n_rejection];
        EEG_epoch_all.dataset          = [EEG_epoch_all.subject, repmat({dataset}, 1, EEG_epoched.trials)];
        


    end % idx_subject = subjects  % iterate over subjects
end % for idx_dataset = 1:length(datasets)

%% Save common struct

% start timer
tic

% start notification
disp('Start saving common struct...')

filepath = [path.save 'prep_1_' type];
% create folder if not existent
if ~exist(filepath, 'dir')
    mkdir(filepath)
end

% Save
save([filepath '\EEG_prep_1_' type '_total'], ...
    "EEG_epoch_all", "n_rejection", '-v7.3');

% timing notification
fprintf('Saving completed. Expired time: %.0f seconds \n', toc)


