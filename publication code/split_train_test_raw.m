%% Settings

clear
clc

path.load = '';
path.save = '';

subjects = 1:50;

datasets = {'young_control', 'stroke_ipsi', 'stroke_contra', 'elderly_control'};

split_time = 3*60;      % training set is first 3 min of rsEEG


%% Initialization

eeglab


%% Iterate over datasets

for idx_dataset = 1:length(datasets)

    dataset = datasets{idx_dataset};
    fprintf('\n \n START DATASET: %s \n \n', dataset);

    for idx_subject = subjects  % iterate over subjects

        formatSpec = '%03.0f';


        %% Load standardized resting state data

        % start timer
        tic

        % start notification
        fprintf('\nStarting participant %03.0f. \n', idx_subject)

        filepath = [path.load '\' dataset '\EEG_rs_compressed_' dataset ...
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


        %% Split data into training and test

        EEG_train   = pop_select(EEG, 'time', [0 split_time]);
        EEG_test    = pop_select(EEG, 'time', [split_time EEG.xmax]);

        clear EEG


        %% Drop participant if too little data is in test set

        if EEG_test.xmax < 60 % if less than one minute of test data (4 min total)
            continue
        end


        %% Save training data

        filepath = [path.save 'raw_train\' dataset];
        % create folder if not existent
        if ~exist(filepath, 'dir')
            mkdir(filepath)
        end

        % Save under variable name "EEG"

        EEG = EEG_train;

        % save 

        save([filepath '\EEG_rs_train_' dataset '_'  num2str(idx_subject,formatSpec)], ...
            "EEG", '-v7.3');

        clear EEG


        %% Save test data

        filepath = [path.save 'raw_test\' dataset];
        % create folder if not existent
        if ~exist(filepath, 'dir')
            mkdir(filepath)
        end

        % Save under variable name "EEG"
        EEG = EEG_test;

        % Save 
        save([filepath '\EEG_rs_test_' dataset '_'  num2str(idx_subject,formatSpec)], ...
            "EEG", '-v7.3');

        clear EEG
        

    end % idx_subject = subjects  % iterate over subjects
end % for idx_dataset = 1:length(datasets)
