function get_ground_truth_phase_estimates(type, settings)
% Load data that has been downsampled, detrended, and laplace transformed
% to estimate the ground truth phase using filtering and a hilbert
% transform.


%% Loading Settings

path.load = '';
path.save = '';


%% Load data preprocessed and data raw (for plotting)

% start timer
tic

% start notification
fprintf('\nStart loading data \n')

filepath = [path.save 'prep_1_' type];
% create folder if not existent
if ~exist(filepath, 'dir')
    mkdir(filepath)
end


% Save
load([filepath '\EEG_prep_1_' type '_total']);

% timing notification
fprintf('loading completed. Expired time: %.0f seconds \n', toc)


%% Filter settings

filtering.order = settings.bp_order;
figure


win_forecast = EEG_epoch_all.t_epoched_raw;    % ms

plotting = true;

filtering.filter = settings.filter;


%% Initialize saving array

phase_all.subject = [];
phase_all.dataset = [];
phase_all.ground_truth = [];

phase_full_all.subject = [];
phase_full_all.dataset = [];
phase_full_all.ground_truth = [];
phase_full_all.time = win_forecast;

EEG_epoch_all.data_cont_filt = {};


%% Iterate over subjects for continuous phase analysis
subject_names = unique(EEG_epoch_all.subject);

for idx_subject = 1:numel(subject_names)

    % retrieve current data
    subject_name_current = subject_names(idx_subject);
    data                 = EEG_epoch_all.data_cont{idx_subject};
    epochs_current       = EEG_epoch_all.events_cont{idx_subject};
    dataset_current      = EEG_epoch_all.dataset(idx_subject);
    time_current         = EEG_epoch_all.t_cont(idx_subject);


    %% Filter

    data = filtfilt(filtering.filter, 1, data);


    %% Phase estimation

    % hilbert transform all data
    hilb = hilbert(data);

    % min(time_current)

    % correct for EEGlab epoching until one sample prior to 0
    % TODO plus one?
    epochs_current = round(epochs_current);
    phase.ground_truth = angle(hilb(epochs_current-1));

    % make a matrix with all locations for the full gt
    locs_epochs = repmat(epochs_current', 1, numel(win_forecast));
    locs_epochs = locs_epochs + repmat(win_forecast, numel(epochs_current), 1);
    phase_full.ground_truth = angle(hilb(locs_epochs));


    %% Save into common array

    % save phase at t = 0
    phase_all.subject = [phase_all.subject, repmat(subject_name_current, 1, length(phase.ground_truth))];
    phase_all.dataset = [phase_all.dataset, repmat({dataset_current}, 1, length(phase.ground_truth))];
    phase_all.ground_truth = [phase_all.ground_truth, phase.ground_truth];
    phase_all.time = 0;

    % save phase at whole epoch
    phase_full_all.subject = [phase_full_all.subject, repmat(subject_name_current, 1, length(phase.ground_truth))];
    phase_full_all.dataset = [phase_full_all.dataset, repmat({dataset_current}, 1, size(phase_full.ground_truth,1))];
    phase_full_all.ground_truth = [phase_full_all.ground_truth; phase_full.ground_truth];
    
    % save filtered data
    EEG_epoch_all.data_cont_filt = [EEG_epoch_all.data_cont_filt {data}];


end


%% Plot raw data

nexttile
hold on; grid on;

percentage = 0.25;

% color alpha
alpha = 0.1;

% Calculate how wide the array of trials is
width_half = 2*pi*percentage/2;
time = EEG_epoch_all.t_epoched_raw;

phase_data = phase_all.ground_truth;

% Select peak data
plot_locs = phase_data < - pi + width_half | phase_data > pi - width_half;
plot_data_all = zscore(EEG_epoch_all.data_epoched_raw, 0, 1);
plot_data = plot_data_all(:,plot_locs);

% plot mu and 95% CI
[mu, CI_lower, CI_upper] = calculate_95_CI(plot_data');
plot(time, mu, 'Color', 'c', 'LineWidth',2, 'DisplayName', 'Trough')
fill([time, fliplr(time)], [CI_upper, fliplr(CI_lower)], ...
    'c', 'EdgeColor', 'none', 'FaceAlpha', alpha);  % light pink fill

% Select trough phase_data
plot_locs = phase_data > -width_half & phase_data < width_half;
plot_data = plot_data_all(:,plot_locs);

% plot 95% CI
[mu, CI_lower, CI_upper] = calculate_95_CI(plot_data');
plot(time, mu, 'Color', 'm', 'LineWidth',2, 'DisplayName', 'Peak')
fill([time, fliplr(time)], [CI_upper, fliplr(CI_lower)], ...
    'm', 'EdgeColor', 'none', 'FaceAlpha', alpha);  % light pink fill

% Select falling phase_data
plot_locs = phase_data > pi/2 - width_half & phase_data < pi/2 + width_half;
plot_data = plot_data_all(:,plot_locs);

% plot 95% CI
[mu, CI_lower, CI_upper] = calculate_95_CI(plot_data');
plot(time, mu, 'Color', [0.75 0.25 1], 'LineWidth',2, 'DisplayName', 'Rising')
fill([time, fliplr(time)], [CI_upper, fliplr(CI_lower)], ...
    [0.75 0.25 1], 'EdgeColor', 'none', 'FaceAlpha', alpha);  % light pink fill

% Select rising phase_data
plot_locs = phase_data > -pi/2 -width_half & phase_data < -pi/2 + width_half;
plot_data = plot_data_all(:,plot_locs);

% plot 95% CI
[mu, CI_lower, CI_upper] = calculate_95_CI(plot_data');
plot(time, mu, 'Color', [0.25 0.75 1], 'LineWidth',2, 'DisplayName', 'Falling')
fill([time, fliplr(time)], [CI_upper, fliplr(CI_lower)], ...
    [0.25 0.75 1], 'EdgeColor', 'none', 'FaceAlpha', alpha);  % light pink fill

% plot specs
ylim([-.5 .5]); xlim([-300 100])
xline(0); yline(0); 

xlabel('time [ms]'); ylabel('amplitude (z-scored) [V]')

l = legend(['Trough' "" "Peak" "" "Falling" "" "Rising" "" ""], 'Location','northoutside');
l.NumColumns = 4;
title("Show binning on ground truth raw data, order = " + num2str(filtering.order));

drawnow

%% Plot distribution of phase bins

nexttile; hold on; grid on

numBins = 13;

% Define bin edges from combined range
edges = linspace(-pi, pi, numBins+1);

% Plot histogram with uniform PDF overlay
histogram(phase_data, edges, 'FaceColor','c', 'FaceAlpha',0.5, 'EdgeAlpha', 0.5)

hold on
yline(length(phase_data)/numBins, 'k', 'LineWidth',2)
xlabel('Data Value');
ylabel('PDF');
xlim([-pi pi])
legend('Ground truth', 'Empirical PDF', 'Uniform PDF', 'Location', 'south');

drawnow

clear h 


%% Save

% start timer
tic

% start notification
disp('Start saving common struct...')

filepath = [path.save 'ground_truth_' type];
% create folder if not existent
if ~exist(filepath, 'dir')
    mkdir(filepath)
end

% Save
save([filepath '\EEG_phase_ground_truth_' type '_total'], ...
    "EEG_epoch_all", "n_rejection", "phase_full_all", "phase_all", '-v7.3');

% timing notification
fprintf('Saving completed. Expired time: %.0f seconds \n', toc)