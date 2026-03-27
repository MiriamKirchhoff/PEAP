function ETP_parameters_optimal = train_ETP(settings)


%% Train ETP phase estimation based on resting state

path.load = '';
path.save = '';

type = 'train';


%% Load ground truth and epoched data

filepath = [path.save 'ground_truth_' type];
load([filepath '\EEG_phase_ground_truth_' type '_total']);

ETP_parameters_optimal = [];

%% Settings

subjects = unique(EEG_epoch_all.subject);

for idx_subject = 1:length(subjects)  % iterate over subjects

    subject_current = subjects(idx_subject);
    disp(['start ' num2str(idx_subject) ' of ' num2str(length(subjects))])
    data = EEG_epoch_all.data_cont_filt{idx_subject};
    fs = EEG_epoch_all.fs;


    %% Estimate phase using ETP

    do_preprocess = false;

    if iscolumn(data)
        data = data';
    end

    [estimated_phase,fullCycle] = ETP_AutoCorrect_edge(data, do_preprocess, settings.bp_frequencies, fs);

    nexttile;
    histogram(estimated_phase,37, "BinLimits",[-pi pi], 'FaceColor','k', 'FaceAlpha',0.3);
    xline(0, 'LineWidth',2)
    drawnow


    %% Save data in array

    ETP_parameters_optimal.estimated_phase{idx_subject} = estimated_phase;
    ETP_parameters_optimal.fullCycle(idx_subject) = fullCycle;
    ETP_parameters_optimal.subject(idx_subject) = subject_current;

end % idx_subject = 1:length(subjects)  % iterate over subjects


%% Save data

filepath = [path.save '\algorithm_training_results\'];

if ~exist(filepath, 'dir')
    mkdir(filepath)
end

save([filepath 'ETP'], 'ETP_parameters_optimal', '-v7.3')

