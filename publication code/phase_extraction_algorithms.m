% Script to get the phase estimates of the preprocessed data

% Script to filter data
%% Settings

path.load = '';
path.save = '';

type = 'test';
target_time = -1;     % [ms]


%% Load ground truth and epoched data

filepath = [path.save 'filtered_' type];
load([filepath '\filtered_' type '_total'])


%% Load individualized settings for phastimate and ETP

filepath = [path.save 'algorithm_training_results\'];
load([filepath 'ETP.mat'])
load([filepath 'phastimate.mat'])
load('settings.mat')
win_forecast = EEG_epoch_all.t_epoched_raw;    % [ms]

%%%%%%%%%%%%%%%%%%%%%%%% For subjects without individualized settings
phast.ord = 30;
phast.hilbertwindow = 128;
phast.offset_correction = 0;
phast.edge = 65;
etp.edge = 40;
%%%%%%%%%%%%%%%%%%%%%%%%


%% Iterate over different types of preprocessing pipelines

% find all fields with filtered data, but not time vectors
temp_fieldnames = fieldnames(EEG_epoch_all);
temp_locations = contains(fieldnames(EEG_epoch_all), 'data_filtered') & ...
    ~contains(fieldnames(EEG_epoch_all), 't_data_filtered');
filter_approach_fieldnames = temp_fieldnames(temp_locations);

clear temp_locations temp_fieldnames

for idx_approach = 1:length(filter_approach_fieldnames)

    approach_fieldname = filter_approach_fieldnames{idx_approach};
    temp = strsplit(approach_fieldname, '_');
    approach_name = temp{end};
    clear temp

    fprintf('Start analysis of %s \n', approach_name)

    % Get data
    data = EEG_epoch_all.(approach_fieldname);
    time = EEG_epoch_all.(['t_' approach_fieldname]);

    % determine target timepoint
    time_target_sample = find(round(time) == target_time);


    %% Phase extraction 1: Hilbert

    disp('start hilbert')

    % exclude NaN trials
    nantrials = logical(sum(~isfinite(data)));
    data_extended_PAR(:, nantrials) = 0;

    hilb = hilbert(data);

    % insert NaNs again
    hilb(:,nantrials) = NaN;
    phase_all.([approach_name 'hilbert']) = angle(hilb(time_target_sample, :));
    phase_full_all.([approach_name 'hilbert']) = angle(hilb)';
    phase_full_all.(['t_' approach_name 'hilbert']) = time;



    %% Phase extraction 2: Phastimate

    disp('start phastimate')

    n_interpolation_points = 100;

    data_phast = data(1:time_target_sample,:);
    phase_full_all.(['t_' approach_name 'phastimate']) = ...
        [time(1:time_target_sample) 0:n_interpolation_points-1];
    phase_full_all.([approach_name 'phastimate']) = ...
        NaN(length(phase_full_all.subject), ...
        length(phase_full_all.(['t_' approach_name 'phastimate'])));

    subject = unique(EEG_epoch_all.subject);
    for idx_subject = 1:length(subject)
        subject_current = subject(idx_subject);
        data_ind = data_phast(:, EEG_epoch_all.subject == subject_current);
        if settings.bp_frequencies(1) == 9
            phast_ind = phastimate_parameters_optimal.settings{ ...
                cell2mat(phastimate_parameters_optimal.subject) == subject_current};

        else % use default settings

            phast_ind = phast;
            phast_ind.optim_ar_ord = phast.ord;
            phast_ind.optim_window = phast.hilbertwindow;
            phast_ind.optim_edge = phast.edge;
        end

        phase_all.([approach_name 'phastimate'])(...
            EEG_epoch_all.subject == subject_current) = phastimate_adapted(...
            data_ind, phast_ind.optim_edge, phast_ind.optim_ar_ord, ...
            phast_ind.optim_window, phast.offset_correction);
        [~, ~, temp] = phastimate_adapted(...
            data_ind, phast_ind.optim_edge, phast_ind.optim_ar_ord, ...
            phast_ind.optim_window, phast.offset_correction, 100 + phast_ind.optim_edge);

        temp = [nan(size(temp, 2), phast_ind.optim_edge) temp'];
        phase_full_all.([approach_name 'phastimate'])(...
            EEG_epoch_all.subject == subject_current, :) = temp;
    end % for idx_subject = 1:length(subject)

    clear temp


    %% Phase extraction 3: ETP

    disp('start etp')

    data_etp = data(1:time_target_sample,:);
    time_etp = time(1:time_target_sample);
    
    min_dist = EEG_epoch_all.fs/(settings.bp_frequencies(2)+1);

    % determine parameters for all subjects and trials
    step_of_phase_per_sample = ...
        2*pi./ETP_parameters_optimal.fullCycle;

    for idx_trial = 1:size(data,2)

        % find subject
        sub = EEG_epoch_all.subject(idx_trial);
        loc = ETP_parameters_optimal.subject == sub;

        step_of_phase_per_sample_current = step_of_phase_per_sample(loc);

        chunk = data_etp(:,idx_trial)';

        locs_hi = mypeakseek(chunk(1:end-etp.edge), min_dist);
        

        % interpolate phase of current cycle
        if ~isempty(locs_hi)

            loc_hi = locs_hi(end);

            % get oscillation
            % cos((time_etp - time_etp(loc_hi)) *  step_of_phase_per_sample_current )
            temp = (win_forecast - time_etp(loc_hi)) *  step_of_phase_per_sample_current;
            phase_full(idx_trial, :) = mod(temp + pi, 2*pi) - pi;
            phase(idx_trial) = phase_full(idx_trial,win_forecast == 0);

        else
            phase(idx_trial)  = NaN;
            phase_full(idx_trial, :) = NaN(size(phase_full(1, :)));
        end


    end % for idx_trial = 1:EEG_epochs.trials

    % correct all to be in [-pi pi]
    phase_all.([approach_name 'etp']) = phase;
    phase_full_all.([approach_name 'etp']) = phase_full;

    phase_full_all.(['t_' approach_name 'etp']) = ...
        win_forecast;

    clear chunk locs_hi t phase phase_full sub loc ...
        step_of_phase_per_sample_current


end % for idx_approach = 1:length(filter_approach_fieldnames)


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
