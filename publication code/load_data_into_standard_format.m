%% Load data and save in standard EEGlab format
% Only keeps data that is relevant for subsequent analysis


%% Settings

path.load_raw = '';
path.save = '';

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


        %% Load raw resting state REFTEP data

        % start timer
        tic

        % start notification
        fprintf('\nStarting participant %03.0f. \n', idx_subject)

        filepath = [dataset '\' num2str(idx_subject,formatSpec)];

        
        % check if datafile exists
        if ~exist(filepath, 'file')
            fprintf('\nNo data found for participant %03.0f. \n', idx_subject)
            continue
        end

        % load EEG data
        disp('Start loading...')
        load(filepath);


        % settings for montage
        switch dataset
            case 'stroke_ipsi'
                % Ipsilaterial hemisphere, as noted in
                % \\storage.neurologie.uni-tuebingen.de\bbnp_lab\Experimental
                % Data\2021-11 REFSTROKE\REFSTROKE_Sessions.xlsx
                ipsi_loc = {'C4' 'C4' 'C4' 'C4' 'C3' 'C3' 'C3' 'C3' '' 'C4'};
                hjorth.channel = ipsi_loc{idx_subject};
                clear ipsi_loc
            case 'stroke_contra'
                % Contralateral hemisphere, as noted in
                % \\storage.neurologie.uni-tuebingen.de\bbnp_lab\Experimental
                % Data\2021-11 REFSTROKE\REFSTROKE_Sessions.xlsx
                contra_loc = {'C3' 'C3' 'C3' 'C3' 'C4' 'C4' 'C4' 'C4' '' 'C3'};
                hjorth.channel = contra_loc{idx_subject};
                clear contra_loc
            otherwise
                hjorth.channel = 'C3';
        end

        % timing notification
        fprintf('loading completed. Expired time: %.0f seconds \n', toc)

        EEG = eeg_checkset(EEG, 'eventconsistency');
        EEG = eeg_checkset(EEG, 'makeur');

        %% Sub-select relevant channels to save processing time

        switch hjorth.channel  % get surrounding variables for hjorth filter
            case 'C3'
                hjorth.electrodes = {'FCC5h', 'FCC3h', 'CCP5h', 'CCP3h'};
            case 'C4'
                hjorth.electrodes = {'FCC6h', 'FCC4h', 'CCP6h', 'CCP4h'};
        end
        temp = struct2cell(EEG.chanlocs);
        hjorth.channel_idx = find(strcmpi(hjorth.channel, temp(1,:)));
        for i = 1:length(hjorth.electrodes)
            hjorth.electrode_idx(i) = find(strcmpi(hjorth.electrodes{i}, temp(1,:)));
        end
        clear i

        EEG = pop_select(EEG, 'channel', ...
            [hjorth.channel_idx, hjorth.electrode_idx]);


        %% Save compressed data

        % create folder if not existent
        if ~exist([path.save '\' dataset], 'dir')
            mkdir([path.save '\' dataset])
        end

        save([path.save '\' dataset ...
            '\EEG_rs_compressed_' dataset '_'  num2str(idx_subject,formatSpec)], ...
            "EEG", '-v7.3');


    end % idx_subject = subjects  % iterate over subjects
end % for idx_dataset = 1:length(datasets)