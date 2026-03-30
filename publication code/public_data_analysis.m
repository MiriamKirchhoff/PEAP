%% Analysis of public data
% This analysis contains the plots that can be re-created using the data
% that can be shared publically.

%% Settings

path.load = '';
path.save = '';

type = 'test';


%% Colors

colors.magenta  = [220 0 90]/255;
colors.blue     = [36 81 219]/255;
colors.cyan     = [0 220 220]/255;
colors.green    = [160 220 90]/255;
colors.orange   = [219 116 0]/255;
colors.purple   = [191 93 219]/255;
colors.grey     = [1 1 1]*0.5;

% set color alpha
f_alpha = 0.3;
colors.f_alpha = f_alpha;
e_alpha = 0.5;


%% Load ground truth and epoched data

filepath = [path.load 'phase_' type];
load('public_data.mat')
% Get fieldnames
approach_fieldnames             = fieldnames(phase_all);
approach_fieldnames([1 2 4 11]) = [];
approach_fieldnames_no_gt       = approach_fieldnames(2:end);
n_approaches                    = length(approach_fieldnames_no_gt);
n_samples                       = length(phase_all.ground_truth);

% Set dataset names, subject labels, subject - dataset interaction
phase_all.subject_dataset       = phase_all.subject;
phase_all.dataset               = floor(phase_all.subject/100);
dataset_names = {'control', 'control', 'Stroke_ipsi', 'Stroke_contra', 'Stroke_control'};
phase_all.dataset_label         = dataset_names(phase_all.dataset);
n_datasets                      = length(unique(phase_all.dataset_label));
dataset_names                   = unique(dataset_names);

% set subjects from stroke data ipsi and contra to same name
loc_stroke = matches(phase_all.dataset_label, 'Stroke_ipsi') | ...
    matches(phase_all.dataset_label, 'Stroke_contra');
phase_all.subject(loc_stroke)   = mod(phase_all.subject(loc_stroke), 100);

phase_full_all.subject_dataset  = phase_all.subject_dataset;
phase_full_all.dataset          = phase_all.dataset;
phase_full_all.dataset_label    = phase_all.dataset_label;
phase_full_all.subject          = phase_all.subject;

subject_dataset_labels          = unique(phase_all.subject_dataset);
n_subject_dataset               = length(subject_dataset_labels);

clear loc_stroke



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                       %%
%%                       COMPUTE MEASURES                                %%
%%                                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create table with all data in long format

indices_nan_errors = zeros(size(phase_all.ground_truth));
% % Find invalid data points
for idx_approach = 1:n_approaches
    indices_nan_errors = indices_nan_errors + isnan(phase_all.(approach_fieldnames_no_gt{idx_approach}));
end
indices_nan_errors = indices_nan_errors > 0;

% Transform data to long format
% Initialize arrays
subject_all = repmat(phase_all.subject(~indices_nan_errors)', n_approaches, 1);
dataset_all = repmat(phase_all.dataset_label(~indices_nan_errors)', n_approaches, 1);
gt_all = repmat(phase_all.ground_truth(~indices_nan_errors)', n_approaches, 1);
approaches_all = repelem(approach_fieldnames_no_gt, sum(~indices_nan_errors));
error_all = [];
estimate_all = [];

% Concatenate all approaches values
for idx_approach = 1:n_approaches
    estimate_all = [estimate_all; double(phase_all.(approach_fieldnames_no_gt{idx_approach})(~indices_nan_errors))'];    
end

data_tbl = table;
data_tbl.subject        = categorical(subject_all);
data_tbl.dataset        = categorical(dataset_all);
data_tbl.approaches     = categorical(approaches_all);
data_tbl.ground_truth   = gt_all;
data_tbl.estimate       = estimate_all;
data_tbl.estimate_sin   = sin(data_tbl.estimate);
data_tbl.estimate_cos   = cos(data_tbl.estimate);


%% Compute errors

data_tbl.error = circ_dist(data_tbl.estimate, data_tbl.ground_truth);
data_tbl.error_normalized = data_tbl.error/pi;
data_tbl.errors_abs     = abs(data_tbl.error);
data_tbl.errors_sqrt    = sqrt(data_tbl.errors_abs);
data_tbl.accuracy_sqrt  = 1-data_tbl.errors_sqrt/sqrt(pi);
data_tbl.accuracy       = (1-data_tbl.errors_abs/pi)*100;


%% Compute SNR
% 
% clear snr_mu
% 
% fs = 1000;
% for idx_trial = 1:size(EEG_epoch_all.data_epoched_raw, 2)
%     % snr_trial(idx_trial) = snr(EEG_epoch_all.data_epoched_raw(:,idx_trial), 1000);
%     % snr_filtered_trial(idx_trial) = snr(EEG_epoch_all.data_filtered_default(:,idx_trial), 1000);
%     [pxx,f] = pwelch(EEG_epoch_all.data_epoched_raw(:,idx_trial),[],[],[],fs);
%     snr_mu(idx_trial) = 10*log10(bandpower(pxx,f,[9 13],'psd') / ...
%         (bandpower(pxx,f,[0 fs/2],'psd') - bandpower(pxx,f,[9 13],'psd')) );
% end
% 
% % data_tbl.snr = repmat(snr_trial(~indices_nan_errors)', n_approaches, 1);
% % data_tbl.snr_filterd = repmat(snr_filtered_trial(~indices_nan_errors)', n_approaches, 1);
% data_tbl.snr_mu = repmat(snr_mu(~indices_nan_errors)', n_approaches, 1);


%% Define bins for analyses

bins.n      = 8;
bins.edges  = linspace(-pi, pi, bins.n+1);
bins.labels = bins.edges(1:end-1);
bins.names = ["trough" "early rising" "rising" "late rising" ...
    "peak" "early falling" "falling" "late falling"];


%% Compute bins for phase estimates

data_tbl.estimate_bin = categorical(get_trough_centered_bin_incides(data_tbl.estimate, bins));

clear temp


%% Toss hilbert transform
% There is no reason to investigate this method as it is not suitable for
% the issue at hand

data_tbl(data_tbl.approaches == "defaulthilbert", :) = [];
data_tbl.approaches = removecats(data_tbl.approaches);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                       %%
%%                          STATISTICS                                   %%
%%                                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MAIN ANALYSIS: Compute linear mixed effects model

data_tbl.subject_dataset = categorical(strcat(string(data_tbl.subject), "_", string(data_tbl.dataset)));

% fit linear mixed effects model:
% fixed effects: estimator, dataset, phase bin
% random effect: subject
% Do not inlcude hilbert due to not all phase bins represented
lme = fitlme(data_tbl, 'accuracy_sqrt ~ approaches * estimate_bin + approaches * dataset + estimate_bin  + (1|subject) + (1|subject_dataset)');

% calculate anova of this model
anovaTbl = anova(lme);
disp(anovaTbl)

clear idx_dataset idx_approach idx_dataset n


%% Get all medians +- mad at t = -1 ms

median_approaches = varfun(@median,data_tbl,'GroupingVariables','approaches','InputVariables',{'accuracy', 'error_normalized'});
mad_approaches = varfun(@mad,data_tbl,'GroupingVariables','approaches','InputVariables',{'accuracy', 'error_normalized'});

group_statistics = join(median_approaches, mad_approaches);
group_statistics.approaches = renamecats(group_statistics.approaches,["PEAPhilbert", "PARhilbert" "PEAPetp" "defaultphastimate" "defaultetp"],["PEAP", "PhastPadding" "PEAP + ETP" "Phastimate" "ETP"]);

group_statistics = sortrows(group_statistics, "median_accuracy", "descend");
disp(group_statistics)
writetable(group_statistics,  [path.save 'accuracies.csv'])

%% Get all contrasts

coefNames = lme.CoefficientNames;

% Get baseline approaches (the reference level automatically chosen by MATLAB)
% MATLAB uses alphabetical ordering, so the first in categories(data_tbl.approaches)
baseline = 'ground_truth';

% Store results
contrast_results = [];

fprintf('Comparing each pair of approaches:\n');

% First, compare each padding vs not padded
padding_options = {'PEAP' 'PAR'};
method_options = {'etp' 'phastimate'};

for idx_method = 1:length(method_options)
    for idx_padding = 1:length(padding_options)
        default = ['default' method_options{idx_method}];
        padded = [padding_options{idx_padding} method_options{idx_method}];

        % Initialize contrast vector
        C = zeros(1, numel(coefNames));

        idx1 = strcmp(coefNames, ['approaches_' default]);
        C(idx1) = 1;
        idx2 = strcmp(coefNames, ['approaches_' padded]);
        C(idx2) = -1;

        % Run contrast test
        [p, F, df1, df2] = coefTest(lme, C);
        contrast_results = [contrast_results; {default, padded, p, F, df1, df2}];
    end
end

% Then, compare list of resulting subset
subset_approaches = approach_fieldnames_no_gt([1 4 3 9 7 8]);
n_subset_approaches = length(subset_approaches);

for idx_approach_1 = 1:n_subset_approaches
    for idx_approach_2 = idx_approach_1+1:n_subset_approaches
        m1 = subset_approaches{idx_approach_1};
        m2 = subset_approaches{idx_approach_2};

        if matches(m1, 'PEAPetp') && matches(m2, 'defaultetp')
            continue; end

        % Initialize contrast vector (same length as number of fixed effects)
        C = zeros(1, numel(coefNames));

        % Get positions of approaches coefficients
        % method_A means difference from baseline, so comparisons must adjust for baseline

        if ~strcmp(m1, baseline)
            idx1 = strcmp(coefNames, ['approaches_' m1]);
            C(idx1) = 1;
        end
        if ~strcmp(m2, baseline)
            idx2 = strcmp(coefNames, ['approaches_' m2]);
            C(idx2) = -1;
        end

        % Run contrast test
        [p, F, df1, df2] = coefTest(lme, C);
        contrast_results = [contrast_results; {m1, m2, p, F, df1, df2}];
    end
end

% Convert to table
contrast_table = cell2table(contrast_results, ...
    'VariableNames', {'Method1', 'Method2', 'pValue', 'Fstat', 'DF1', 'DF2'});

[~, ~, ~, contrast_table.pValue_corrected] = fdr_bh(contrast_table.pValue);

% Round values
contrast_table.pValue = round(contrast_table.pValue, 3);
contrast_table.Fstat = round(contrast_table.Fstat, 3);
contrast_table.pValue_corrected = round(contrast_table.pValue_corrected, 3);
disp(contrast_table)

% Extract which approaches were improved by preprocessing pipelines -> Can
% be used for plotting, only plot if improved
approaches_processing = ["default", "PAR", "PEAP"];
approaches_phase_extraction = ["hilbert", "etp", "phastimate"];
approaches_significant = approach_fieldnames_no_gt;

for idx_approach_phase = 1:length(approaches_phase_extraction)
    for idx_approach_processing = 2:length(approaches_processing)

        % exclude if has been compared before
        
        
        % default name
        def_name = 'default' + ...
            approaches_phase_extraction(idx_approach_phase);

        % alternative name
        alt_name = approaches_processing(idx_approach_processing) + ...
            approaches_phase_extraction(idx_approach_phase);

        % find default positions
        def_pos = contains(contrast_table.Method1, def_name) | ...
            contains(contrast_table.Method2, def_name);

        % find alternative position
        alt_pos = contains(contrast_table.Method1, alt_name) | ...
            contains(contrast_table.Method2, alt_name);

        % find if contrast is significant
        pos = def_pos & alt_pos;

        % if not significantly improved, remove from list
        if contrast_table.pValue_corrected(pos) > 0.05
            approaches_significant( ...
                matches(approaches_significant, alt_name)) = [];
        end
    end
end

clear p F df1 df2 m1 m2 p F df1 df2 C idx2 idx1 idx_dataset j  baseline anovaTbl...
    def_pos def_name alt_name alt_pos pos

% Export to excel
writetable(contrast_table(:, [1 2 7 4:6]), [path.save 'multiple _comparisons.csv']);
disp('multiple comparisons table successfully exported')


%% Prepare data for figures

idx_plot = contains(string(data_tbl.approaches), approaches_significant);
data_tbl_plot = data_tbl(idx_plot,:);
% remove extra categories
data_tbl_plot.approaches = categorical(string(data_tbl_plot.approaches));

% Sort trials

data_tbl_plot.dataset = renamecats(data_tbl_plot.dataset,["control", "Stroke_control" "Stroke_contra" "Stroke_ipsi"],["Young control", "Elderly control" "Stroke intact" "Stroke affected"]);
data_tbl_plot.dataset = reordercats(data_tbl_plot.dataset, ["Young control", "Elderly control" "Stroke intact" "Stroke affected"]);

%data_tbl_plot.approaches = renamecats(data_tbl_plot.approaches,["PEAPhilbert", "PARhilbert" "defaultphastimate" "defaultetp"],["PEAP", "PhastPadding" "Phastimate" "ETP"]);
%data_tbl_plot.approaches = reordercats(data_tbl_plot.approaches, ["PEAP", "PhastPadding" "PEAP + Phastimate" "Phastimate", "SSPE" "ETP"]);

data_tbl_plot.approaches = renamecats(data_tbl_plot.approaches,["PEAPhilbert", "PARhilbert" "PEAPetp" "defaultphastimate" "defaultetp"],["PEAP", "PhastPadding" "PEAP + ETP" "Phastimate" "ETP"]);
data_tbl_plot.approaches = reordercats(data_tbl_plot.approaches, ["PEAP" "PhastPadding" "PEAP + ETP"  "SSPE" "Phastimate" "ETP"]);


%% PAPER FIGURE 3a: Accuracy

figure
t = tiledlayout('flow', "TileSpacing", "compact", 'Padding', 'tight');
title(t, "Accuracy across methods and datasets", 'FontName', 'Times');


nexttile()
b = boxchart(data_tbl_plot.approaches, data_tbl_plot.accuracy, 'GroupByColor', data_tbl_plot.dataset, 'ColorGroupWidth',0.8);

boxplot_settings(b, colors);

legend('Interpreter','none', 'Location','southeast')
ylim([0, 100])
ylabel('Accuracy [%]');
grid on;

l = legend('Interpreter','none', 'Box','off');
l.Location = "northoutside";
l.NumColumns = 2;
l.Box = "on";

clear t b 

fontsize(gcf,scale=3)
set(gca,'FontName', 'Times')
drawnow

%% PAPER FIGURE 3b: Comparison of forecasting

figure
t = tiledlayout("TileSpacing", "compact", 'Padding', 'tight');
title(t, "Forecasting accuracy", 'FontName', 'Times');

nexttile([2,2]);hold on; grid on

time_min = -100;
time_max = 50;
time = time_min:time_max;
gt = phase_full_all.ground_truth(...
    :, phase_full_all.time >= time_min & phase_full_all.time <= time_max);

col = ([...
    colors.magenta; ...
    colors.orange; ...
    colors.green; ...
    colors.cyan; ...
    colors.blue; ...
    colors.purple] ...
    ) ;

idx_c = 1;

for idx_approach = 1:numel(approaches_significant)
    % SSPE does not predict, so skip    
    if matches(approaches_significant{idx_approach}, 'SSPE')
        continue
    end

    % find time indices
    temp_time = phase_full_all.(['t_' approaches_significant{idx_approach}]);

    % Select subset of time indices
    time_loc = temp_time >= time_min & temp_time <= time_max;
    acc_temp = (1-abs(circ_dist(phase_full_all.(approaches_significant{idx_approach})(:, time_loc), gt))/pi)*100;
    time_selected = time;

    [mu, CI_lower, CI_upper] = calculate_95_CI_median(acc_temp, 100);
    
    plot(time_selected, mu, 'Color', col(idx_c,:), 'LineWidth',2, 'DisplayName', get_fig_label_approach(approaches_significant{idx_approach}))
    fill([time_selected, fliplr(time_selected)], [CI_upper, fliplr(CI_lower)], ...
        col(idx_c,:), 'EdgeColor', 'none', 'FaceAlpha', f_alpha, 'DisplayName',"\pm 95% CI");  % light pink fill
    idx_c = idx_c + 1;

    xlim([time_min time_max])
    ylim([0.65 1]*100)

    set(gca,'FontName', 'Times')
end
% title('Total')
ax = gca;
h = ax.Children;                       % get all plotted objects
h_even = h(1:2:end);                   % select every second
set(h_even, 'HandleVisibility', 'off') % hide from legend
l = legend('AutoUpdate','off', 'Location', 'south');
xlabel("Time [ms]")
ylabel("Accuracy [%]")
xticks(-100:20:50)
yticks(60:5:100)
xline(0, '-k', 'TMS pulse', LineWidth=2, FontName='Times')

l = legend('Interpreter','none', 'Box','off');
l.Location = "northoutside";
l.NumColumns = 3;
l.Box = "on";


fontsize(gcf,scale=3)
drawnow

clear ax h h_even idx_c l col


%% PAPER FIGURE 4 + ANALYSIS: Uniformity of samples

% bar plot

dataset_current = 'total';

fig = figure(Units="normalized", Position=[.3 .3 .7 .6]);
t = tiledlayout(3, numel(approaches_significant), "TileSpacing", "compact", 'Padding', 'tight');
title(t, "Comparing frequency, accuracy, and bias across estimated phases", 'FontName', 'Times');
results.(dataset_current).sample_uniformity = table('Size', [n_approaches, 4], ...
    'VariableTypes', {'string', 'double', 'double', 'logical'}, ...
    'VariableNames', {'Estimator', 'KS_Stat', 'pValue', 'HReject'});

% Loop over approaches
for idx_approach = 1:numel(approaches_significant)
    approach_current = approaches_significant{idx_approach};
    idx_current = (data_tbl.approaches == approach_current);
    data = data_tbl.estimate(idx_current);

    % compare ground truth distribution to estimated distribution
    [pval, k] = circ_kuipertest(data, phase_all.ground_truth, 1000, 0);

    % Store results
    results.(dataset_current).sample_uniformity.Estimator(idx_approach) = approach_current;
    results.(dataset_current).sample_uniformity.KS_Stat(idx_approach) = k;
    results.(dataset_current).sample_uniformity.pValue(idx_approach) = pval;

    nexttile; hold on; grid on % set(gca, 'YGrid', 'on');

    h = plot_bin_histogram(...
        data, phase_all.ground_truth, bins, pval, colors);

    ylim([0 20])

    title(sprintf('%s', get_fig_label_approach(approach_current)), 'Interpreter','none');
    xlim([0.5 8.5]);
    yticks(0:5:20)
    
    if idx_approach == 1;  ylabel('Percentage [%]'); else; yticklabels([]); end
    set(gca,'FontName', 'Times'); 

end

l = legend('Ground truth', 'Observed', 'Uniform', 'Location', 'eastoutside');
%l.Layout.Tile = 'North';
%l.NumColumns = 3;

% Display results table
disp(results.(dataset_current).sample_uniformity)
drawnow

clear k idx_approach idx_dataset dataset_current dataset_in fig pval


% FIGURE: Compare phase bin error sliding window

for idx_approach = 1:numel(approaches_significant)
    approach_current = approaches_significant{idx_approach};
    data_tbl_current = data_tbl(data_tbl.approaches == approach_current, :);

    disp(get_fig_label_approach(approach_current))

    x = -pi:0.1:pi;
    width_window = 2*pi/bins.n;

    % Use sliding window to calculate mean +- CI for each datapoint

    for idx_x = 1:length(x)

        x_current = x(idx_x);
        x_min = x_current - .5*width_window;
        x_max = x_current + .5*width_window;

        % select all relevant datapoints
        if x_min > -pi && x_max < pi
            data_idx = data_tbl_current.estimate >= x_min &...
                data_tbl_current.estimate < x_max;

        elseif x_min < -pi && x_max < pi
            data_idx = data_tbl_current.estimate >= x_min &...
                data_tbl_current.estimate < x_max;
            % Add values at upper end of the spectrum
            data_idx = data_idx | data_tbl_current.estimate >= x_min + 2*pi;

        elseif x_min > -pi && x_max > pi 
            data_idx = data_tbl_current.estimate >= x_min &...
                data_tbl_current.estimate < x_max;
            % Add values at lower end of the spectrum
            data_idx = data_idx | data_tbl_current.estimate < x_max - 2*pi;
        end

        data = data_tbl_current.accuracy(data_idx);
        mean_x(idx_x) = median(data);
        ci_x(:, idx_x) = bootci(1000,@median,data);

        % plot(data_tbl_current.estimate(data_idx), zeros(1, length(data_tbl_current.estimate(data_idx))), '.')
        % xlim([-pi pi])
        % drawnow



    end

    nexttile; grid on; hold on;

    fill([x'; flipud(x')], [ci_x(2,:)'; flipud(ci_x(1,:)')], ...
         colors.green, 'FaceAlpha', f_alpha, 'EdgeColor', colors.green, 'EdgeAlpha', 1)
    plot(x, mean_x, 'k', 'LineWidth', 2)
    % title(get_fig_label_approach(approach_current))
    xlim([-pi pi])
    ylim([0.69 0.87]*100)
    set(gca,'FontName', 'Times')

    % plot median accuracy across phases

    yline(median(data_tbl_current.accuracy), ':k', 'LineWidth', 2)

    if idx_approach == 1; ylabel('Accuracy [%]');
    else yticklabels([]); end

    xticks(bins.edges(1:2:7))
    xticklabels(["" "" "" ""])

    yticks((0.7:0.05:0.9)*100)
    drawnow
end

l = legend('95% CI', 'Window Median', 'Total Median', 'Location', 'eastoutside');
l.Direction = 'reverse';

clear x mean_x ci_x


for idx_approach = 1:numel(approaches_significant)
    approach_current = approaches_significant{idx_approach};
    data_tbl_current = data_tbl(data_tbl.approaches == approach_current, :);

    disp(get_fig_label_approach(approach_current))

    x = -pi:0.1:pi;
    width_window = 2*pi/bins.n;

    % Use sliding window to calculate mean +- CI for each datapoint

    for idx_x = 1:length(x)

        x_current = x(idx_x);
        x_min = x_current - .5*width_window;
        x_max = x_current + .5*width_window;

        % select all relevant datapoints
        if x_min > -pi && x_max < pi
            data_idx = data_tbl_current.estimate >= x_min &...
                data_tbl_current.estimate < x_max;

        elseif x_min < -pi && x_max < pi
            data_idx = data_tbl_current.estimate >= x_min &...
                data_tbl_current.estimate < x_max;
            % Add values at upper end of the spectrum
            data_idx = data_idx | data_tbl_current.estimate >= x_min + 2*pi;

        elseif x_min > -pi && x_max > pi 
            data_idx = data_tbl_current.estimate >= x_min &...
                data_tbl_current.estimate < x_max;
            % Add values at lower end of the spectrum
            data_idx = data_idx | data_tbl_current.estimate < x_max - 2*pi;
        end

        % Get normalized data in percent
        data = data_tbl_current.error(data_idx)/pi*100;
        mean_x(idx_x) = median(data);
        ci_x(:, idx_x) = bootci(1000,@median,data);

        % plot(data_tbl_current.estimate(data_idx), zeros(1, length(data_tbl_current.estimate(data_idx))), '.')
        % xlim([-pi pi])
        % drawnow

    end

    nexttile; grid on; hold on;

    fill([x'; flipud(x')], [ci_x(2,:)'; flipud(ci_x(1,:)')], ...
         colors.magenta, 'FaceAlpha', f_alpha, 'EdgeColor', colors.magenta, 'EdgeAlpha', e_alpha)
    plot(x, mean_x, 'k', 'LineWidth', 2)
    % title(get_fig_label_approach(approach_current))
    xlim([-pi pi])
    ylim([-15 15])
    yline(0, ':k', 'linewidth', 2)
        set(gca,'FontName', 'Times')

   if idx_approach == 1; ylabel('Bias [%]'); 
   else yticklabels([]); end

   xticks(bins.edges(1:2:7))
   xticklabels(bins.names(1:2:7))
   xtickangle(90)
   drawnow
end
l = legend('95% CI', 'Window Median', 'Location', 'eastoutside');
l.Direction = 'reverse';

fontsize(gcf,scale=2.8); drawnow

clear x mean_x ci_x


%% PAPER FIGURE 3c: Raw data

% The data required for this figure can not be provided publically. For
% more information, please contact Prof. Ziemann.

% figure
% t = tiledlayout(5,10, "TileSpacing", "compact", 'Padding', 'tight');
% title(t, "Binning on unfiltered data", 'FontName', 'Times');
% 
% percentage = 0.25;
% % Calculate how wide the array of trials is
% width_half = 2*pi*percentage/2;
% time = EEG_epoch_all.t_epoched_raw;
% approaches_temp = [{'ground_truth'}; approaches_significant];
% 
% dataset_indices = ones(size(phase_all.dataset));
% for idx_approach = 1:length(approaches_temp)
%     if idx_approach < 3
%         nexttile([3,5]);
%     else
%         nexttile([2,2])
%     end
%     grid on; hold on
% 
%     data = phase_all.(approaches_temp{idx_approach});
%     x = time;
% 
% 
%     % Select trough data
%     plot_locs = data < - pi + width_half | data > pi - width_half;
%     plot_data_all = zscore(EEG_epoch_all.data_epoched_raw, 0, 1);
%     plot_data = plot_data_all(:,plot_locs);
% 
%     % plot mu and 95% CI
%     [mu, CI_lower, CI_upper] = calculate_95_CI(plot_data');
%     plot(time, mu, 'Color', colors.cyan, 'LineWidth',2, 'DisplayName', 'Trough')
% 
%     fill([x, fliplr(x)], [CI_upper, fliplr(CI_lower)], ...
%         colors.cyan, 'EdgeColor', colors.cyan, 'EdgeAlpha', e_alpha, 'FaceAlpha', f_alpha);  % light pink fill
% 
% 
%     % Select rising data
%     plot_locs = data > -pi/2 -width_half & data < -pi/2 + width_half;
%     plot_data = plot_data_all(:,plot_locs);
% 
%     % plot 95% CI
%     [mu, CI_lower, CI_upper] = calculate_95_CI(plot_data');
%     plot(time, mu, 'Color', colors.green, 'LineWidth',2, 'DisplayName', 'Falling')
%     fill([x, fliplr(x)], [CI_upper, fliplr(CI_lower)], ...
%         colors.green, 'EdgeColor', colors.green, 'EdgeAlpha', e_alpha, 'FaceAlpha', f_alpha);  % light pink fill
% 
% 
%     % Select peak data
%     plot_locs = data > -width_half & data < width_half;
%     plot_data = plot_data_all(:,plot_locs);
% 
%     % plot 95% CI
%     [mu, CI_lower, CI_upper] = calculate_95_CI(plot_data');
%     plot(time, mu, 'Color', colors.magenta, 'LineWidth',2, 'DisplayName', 'Peak')
%     fill([x, fliplr(x)], [CI_upper, fliplr(CI_lower)], ...
%         colors.magenta, 'EdgeColor', colors.magenta, 'EdgeAlpha', e_alpha, 'FaceAlpha', f_alpha);  % light pink fill
% 
% 
%     % Select falling data
%     plot_locs = data > pi/2 - width_half & data < pi/2 + width_half;
%     plot_data = plot_data_all(:,plot_locs);
% 
%     % plot 95% CI
%     [mu, CI_lower, CI_upper] = calculate_95_CI(plot_data');
%     plot(time, mu, 'Color', colors.orange, 'LineWidth',2, 'DisplayName', 'Rising')
%     fill([x, fliplr(x)], [CI_upper, fliplr(CI_lower)], ...
%         colors.orange, 'EdgeColor', colors.orange, 'EdgeAlpha', e_alpha, 'FaceAlpha', f_alpha);  % light pink fill
% 
% 
%     % plot specs
%     ylim([-.65 .65]); xlim([-200 65])
%     xline(-1, 'Linewidth', 2); yline(0, 'Linewidth', 2);
%     title(get_fig_label_approach(approaches_temp{idx_approach}));
%     if idx_approach == 1 || idx_approach == 3; ylabel('Amplitude [V]'); 
%     else yticklabels([]); end
%     if idx_approach > 2; yticks(-0.5:0.5:0.5); end
%     if idx_approach == 5; xlabel('Time [ms]'); end
%     set(gca,'FontName', 'Times')
% 
% end
% 
% l = legend(['Trough' "" "Rising" "" "Peak" "" "Falling" "" ""]);
% l.Layout.Tile = 'North';
% l.NumColumns = 4;
% 
% fontsize(gcf,scale=2.8)
% drawnow
% 
% clear mu CI_lower CI_upper l plot_data plot_locs x time width_half t ...
%     percentage


%% Quality control

fig = figure(Units="normalized", Position=[0.5 0.4 0.25 0.25]);
t = tiledlayout('flow', 'TileSpacing','compact');
title(t, 'Power-frequency spectrum', 'FontName', 'timesnewroman')

% testing set
% compute PSD for all subjects


disp('testing quality control')
snr = snr_mu.test;

N  = numel(snr);
mu = mean(snr);
sd = std(snr);
ci = tinv([0.025 0.975], N-1) * sd/sqrt(N);

fprintf('Mean SNR = %.2f dB (95%% CI: [%.2f, %.2f])\n', mu, mu+ci(1), mu+ci(2));

snr = snr_mu.test_filt;

N  = numel(snr);
mu = mean(snr);
sd = std(snr);
ci = tinv([0.025 0.975], N-1) * sd/sqrt(N);
fprintf('Mean filtered SNR = %.2f dB (95%% CI: [%.2f, %.2f])\n', mu, mu+ci(1), mu+ci(2));

% Plot PSD
f_plot   = f(f > 1 & f < 45);
psd_plot = log(psd_all.test((f > 1 & f < 45), :));

nexttile; hold on; grid on


plot(f_plot, mean(psd_plot,2), 'Color', colors.cyan, 'LineWidth', 2, 'DisplayName', 'Test')


% train set

disp('training quality control')
snr = snr_mu.train;

N  = numel(snr);
mu = mean(snr);
sd = std(snr);
ci = tinv([0.025 0.975], N-1) * sd/sqrt(N);

fprintf('Mean SNR = %.2f dB (95%% CI: [%.2f, %.2f])\n', mu, mu+ci(1), mu+ci(2));

snr = snr_mu.train_filt;

N  = numel(snr);
mu = mean(snr);
sd = std(snr);
ci = tinv([0.025 0.975], N-1) * sd/sqrt(N);
fprintf('Mean filtered SNR = %.2f dB (95%% CI: [%.2f, %.2f])\n', mu, mu+ci(1), mu+ci(2));

% Plot PSD
f_plot   = f(f > 1 & f < 45);
psd_plot = log(psd_all.train((f > 1 & f < 45), :));


% plot(f_plot, mean(psd_plot,2), 'Color', 'w', 'LineWidth', 2, 'DisplayName', '\textbf{Training}')
hold on
plot(f_plot, mean(psd_plot,2), 'Color', colors.magenta, 'LineWidth', 2, 'DisplayName', 'Training')


xlim([4 45])
yl = ylim;
patch([9 13 13 9], ...
      [yl(1) yl(1) yl(2) yl(2)], ...
      colors.grey, ...
      'FaceAlpha', 0.15, ...
      'EdgeColor', 'none', ...
      'HandleVisibility','off');

xlabel('Frequency [Hz]')
ylabel('Power [log(dB)]')
legend('NumColumns', 2, 'Interpreter','latex')
fontsize(gcf,scale=1.4)
set(gca,'FontName', 'Times')




%% Functions

function [b_settings] = boxplot_settings(b, colors)

%colors = ['g', 'c', 'b', 'm', 'w'];
if length(b) > 1
    colors_current = {colors.magenta, colors.orange, colors.green, colors.cyan};
else
    colors_current = {colors.grey};
end

for i = 1:length(b)
    b(i).BoxEdgeColor = 'k';
    b(i).BoxFaceColor = colors_current{i};
    b(i).BoxFaceAlpha = colors.f_alpha;
    b(i).MarkerStyle = '.';
    b(i).MarkerColor = 'k';
    b(i).MarkerSize = 1;
    b(i).JitterOutliers = "on";
    b(i).BoxWidth = .8;
    b(i).LineWidth = 2;
end
b_settings = b;

end

function h = plot_bin_histogram(data, ground_truth, bins, pval, colors)

% Add ground truth as stairs
bin_numbers = get_trough_centered_bin_incides(ground_truth, bins);

histogram(bin_numbers, 'DisplayStyle','stairs', 'EdgeColor', 'k', 'EdgeAlpha', 1, 'LineStyle', ':', 'LineWidth',2, 'Normalization','percentage')

bin_numbers = get_trough_centered_bin_incides(data, bins);

% Plot histogram with uniform PDF overlay
h = histogram(bin_numbers, 'FaceColor', colors.cyan, 'FaceAlpha', colors.f_alpha, 'EdgeAlpha', 0.6, 'LineWidth',2, 'Normalization','percentage');

xticks(1:2:bins.n);
names = "";
names(2:2:end) = [];
xticklabels(names)

hold on
yline(100/bins.n, 'k', 'LineWidth',2)
end

function name = get_fig_label_approach(approach)
switch approach
    case 'ground_truth'
        name = 'Ground truth';
    case 'PEAPhilbert'
        name = 'PEAP';
    case 'PEAPphastimate'
        name = 'PEAP + Phastimate';
    case 'PEAPetp'
        name = 'PEAP + ETP';
    case 'PARhilbert'
        name = 'PhastPadding';
    case 'PARphastimate'
        name = 'PhastPadding + Phastimate';
    case 'PARetp'
        name = 'PhastPadding + ETP';
    case 'SSPE'
        name = 'SSPE';
    case 'defaulthilbert'
        name = 'Hilbert';
    case 'defaultphastimate'
        name = 'Phastimate';
    case 'defaultetp'
        name = 'ETP';
end
end

function bin_numbers = get_trough_centered_bin_incides(data, bins)

% center it around trough by shifting the values slightly positively by
% half a bin
temp = data + (2*pi)/bins.n/2;
% move values above pi to trough bin
temp(temp > pi) = -pi + 0.0001;
bin_numbers = discretize(temp, bins.edges);

end

function [fig_new, l] = setings_figures_publications(fig)

% Legend
l = legend('Interpreter', 'none');
l.Layout.Tile = 'North';

% Font
fontsize(gcf,scale=1.4)
set(gca,'FontName', 'Times')

fig_new = fig;

end

