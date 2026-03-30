% Script for analysis 
% Script to get the phase estimates of the preprocessed data
clear all

% Script to filter data
%% Settings

path.load = 'C:\Users\Miriam_Kirchhoff\Documents\MATLAB\REFTEP_preprocessing\Phase estimation scripts\data\';
path.save = 'C:\Users\Miriam_Kirchhoff\Documents\MATLAB\REFTEP_preprocessing\Phase estimation scripts\Results\';


type = 'test';

% dataset_names = {'LOWERLIMB', 'REFTEP', 'REFSTROKE_ipsi', 'REFSTROKE_contra', 'REFSTROKE_healthy'};
dataset_names = {'control', 'control', 'Stroke_ipsi', 'Stroke_contra', 'Stroke_control'};

plot_paper_only = true;

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
load([filepath '\phase_' type '_total'])

% Get fieldnames
approach_fieldnames             = fieldnames(phase_all);
approach_fieldnames([1 2 4 11]) = [];
approach_fieldnames_no_gt       = approach_fieldnames(2:end);
n_approaches                    = length(approach_fieldnames_no_gt);
n_samples                       = length(phase_all.ground_truth);

% Set dataset names, subject labels, subject - dataset interaction
phase_all.subject_dataset       = phase_all.subject;
phase_all.dataset               = floor(phase_all.subject/100);
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

clear snr_mu

fs = 1000;
for idx_trial = 1:size(EEG_epoch_all.data_epoched_raw, 2)
    % snr_trial(idx_trial) = snr(EEG_epoch_all.data_epoched_raw(:,idx_trial), 1000);
    % snr_filtered_trial(idx_trial) = snr(EEG_epoch_all.data_filtered_default(:,idx_trial), 1000);
    [pxx,f] = pwelch(EEG_epoch_all.data_epoched_raw(:,idx_trial),[],[],[],fs);
    snr_mu(idx_trial) = 10*log10(bandpower(pxx,f,[9 13],'psd') / ...
        (bandpower(pxx,f,[0 fs/2],'psd') - bandpower(pxx,f,[9 13],'psd')) );
end

% data_tbl.snr = repmat(snr_trial(~indices_nan_errors)', n_approaches, 1);
% data_tbl.snr_filterd = repmat(snr_filtered_trial(~indices_nan_errors)', n_approaches, 1);
data_tbl.snr_mu = repmat(snr_mu(~indices_nan_errors)', n_approaches, 1);


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


%% PAPER FIGURE 1: Accuracy

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

%% PAPER FIGURE XX: Comparison of forecasting

figure
t = tiledlayout("TileSpacing", "compact", 'Padding', 'tight');
title(t, "Forecasting accuracy", 'FontName', 'Times');

nexttile([2,2]);hold on; grid on

time_min = -100;
time_max = 50;
time = time_min:time_max;
gt = phase_full_all.ground_truth(...
    :, phase_full_all.time >= time_min & phase_full_all.time <= time_max);

% colors = colormap(hsv(numel(approaches_significant)));

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
exportgraphics(fig, 'REFTEP_preprocessing\Phase estimation scripts\figures\bias.tif', 'Resolution', 500);

clear x mean_x ci_x


%% PAPER FIGURE: Raw data

figure
t = tiledlayout(5,10, "TileSpacing", "compact", 'Padding', 'tight');
title(t, "Binning on unfiltered data", 'FontName', 'Times');

percentage = 0.25;
% Calculate how wide the array of trials is
width_half = 2*pi*percentage/2;
time = EEG_epoch_all.t_epoched_raw;
approaches_temp = [{'ground_truth'}; approaches_significant];

dataset_indices = ones(size(phase_all.dataset));
for idx_approach = 1:length(approaches_temp)
    if idx_approach < 3
        nexttile([3,5]);
    else
        nexttile([2,2])
    end
    grid on; hold on

    data = phase_all.(approaches_temp{idx_approach});
    x = time;


    % Select trough data
    plot_locs = data < - pi + width_half | data > pi - width_half;
    plot_data_all = zscore(EEG_epoch_all.data_epoched_raw, 0, 1);
    plot_data = plot_data_all(:,plot_locs);

    % plot mu and 95% CI
    [mu, CI_lower, CI_upper] = calculate_95_CI(plot_data');
    plot(time, mu, 'Color', colors.cyan, 'LineWidth',2, 'DisplayName', 'Trough')
    
    fill([x, fliplr(x)], [CI_upper, fliplr(CI_lower)], ...
        colors.cyan, 'EdgeColor', colors.cyan, 'EdgeAlpha', e_alpha, 'FaceAlpha', f_alpha);  % light pink fill


    % Select rising data
    plot_locs = data > -pi/2 -width_half & data < -pi/2 + width_half;
    plot_data = plot_data_all(:,plot_locs);

    % plot 95% CI
    [mu, CI_lower, CI_upper] = calculate_95_CI(plot_data');
    plot(time, mu, 'Color', colors.green, 'LineWidth',2, 'DisplayName', 'Falling')
    fill([x, fliplr(x)], [CI_upper, fliplr(CI_lower)], ...
        colors.green, 'EdgeColor', colors.green, 'EdgeAlpha', e_alpha, 'FaceAlpha', f_alpha);  % light pink fill


    % Select peak data
    plot_locs = data > -width_half & data < width_half;
    plot_data = plot_data_all(:,plot_locs);

    % plot 95% CI
    [mu, CI_lower, CI_upper] = calculate_95_CI(plot_data');
    plot(time, mu, 'Color', colors.magenta, 'LineWidth',2, 'DisplayName', 'Peak')
    fill([x, fliplr(x)], [CI_upper, fliplr(CI_lower)], ...
        colors.magenta, 'EdgeColor', colors.magenta, 'EdgeAlpha', e_alpha, 'FaceAlpha', f_alpha);  % light pink fill


    % Select falling data
    plot_locs = data > pi/2 - width_half & data < pi/2 + width_half;
    plot_data = plot_data_all(:,plot_locs);

    % plot 95% CI
    [mu, CI_lower, CI_upper] = calculate_95_CI(plot_data');
    plot(time, mu, 'Color', colors.orange, 'LineWidth',2, 'DisplayName', 'Rising')
    fill([x, fliplr(x)], [CI_upper, fliplr(CI_lower)], ...
        colors.orange, 'EdgeColor', colors.orange, 'EdgeAlpha', e_alpha, 'FaceAlpha', f_alpha);  % light pink fill


    % plot specs
    ylim([-.65 .65]); xlim([-200 65])
    xline(-1, 'Linewidth', 2); yline(0, 'Linewidth', 2);
    title(get_fig_label_approach(approaches_temp{idx_approach}));
    if idx_approach == 1 || idx_approach == 3; ylabel('Amplitude [V]'); 
    else yticklabels([]); end
    if idx_approach > 2; yticks(-0.5:0.5:0.5); end
    if idx_approach == 5; xlabel('Time [ms]'); end
    set(gca,'FontName', 'Times')

end

l = legend(['Trough' "" "Rising" "" "Peak" "" "Falling" "" ""]);
l.Layout.Tile = 'North';
l.NumColumns = 4;

fontsize(gcf,scale=2.8)
drawnow

clear mu CI_lower CI_upper l plot_data plot_locs x time width_half t ...
    percentage

%% PLOT: SNR

fig = figure(Units="normalized", Position=[.5 .2 .4 .7]);
t = tiledlayout(2, numel(approaches_significant)/2, "TileSpacing", "compact");
title(t, 'Correlation of SNR with bias and accuracy', 'FontName', 'Times')

% colors_current = {colors.magenta, colors.green, colors.cyan, colors.blue};
colors_current = {"k", "k", "k", "k"};
for idx_approach = 1:length(approaches_significant)
    approach_current = approaches_significant{idx_approach};
    idx_current = (data_tbl.approaches == approach_current);
    snr_current = data_tbl.snr_mu(idx_current);
    bias_current = data_tbl.error(idx_current);
    accuracy = data_tbl.accuracy(idx_current);

    nexttile; grid on; hold on
    set(gca,'FontName', 'Times')

    
    % plot error correlation with fit
    for idx_dataset = 1:n_datasets
        idx_dataset_loc = idx_current & data_tbl.dataset == dataset_names{idx_dataset};
        scatter(data_tbl.snr_mu(idx_dataset_loc), data_tbl.accuracy(idx_dataset_loc), 3, ...
            'MarkerFaceColor', colors_current{idx_dataset}, ...
            'MarkerEdgeColor', colors_current{idx_dataset}, ...
            'MarkerFaceAlpha', 0.6, ...
            'MarkerEdgeAlpha', 0, ...
            'DisplayName', dataset_names{idx_dataset});
    end
    %plot(snr_subject, acc_subject, '.m', 'DisplayName', 'Abs error', 'MarkerSize', 10)
    [correlation_error, pval] = corr(snr_current, accuracy);

    p = polyfit(snr_current, accuracy, 1);
    y_data = polyval(p, snr_current);
    plot(snr_current, y_data, 'k', 'LineWidth', 2)

    if idx_approach == 1; ylabel(t, 'Accuracy [%]', 'FontName', 'Times'); end
    xlabel(t, 'SNR [db]', 'FontName', 'Times')
    %ylim([.45 .95]);
    xlim([-20 6.3])

    title([get_fig_label_approach(approach_current) " \rho = " + num2str(round(correlation_error, 2))])

    fprintf('Correlation SNR and accuracy for %s: rho = %.3f, p = %.3f \n', ...
        get_fig_label_approach(approach_current), ...
        correlation_error, pval)

end

for idx_approach = 1:length(approaches_significant)
    approach_current = approaches_significant{idx_approach};
    idx_current = (data_tbl.approaches == approach_current);
    snr_current = data_tbl.snr_mu(idx_current);
    bias_current = data_tbl.error(idx_current)/pi;
    accuracy = data_tbl.accuracy(idx_current);

    % nexttile; grid on; hold on
    % set(gca,'FontName', 'Times')
    % 
    % 
    % % plot error correlation with fit
    % for idx_dataset = 1:n_datasets
    %     idx_dataset_loc = idx_current & data_tbl.dataset == dataset_names{idx_dataset};
    %     scatter(data_tbl.snr(idx_dataset_loc), data_tbl.error(idx_dataset_loc)/pi, 3, ...
    %         'MarkerFaceColor', colors_current{idx_dataset}, ...
    %         'MarkerEdgeColor', colors_current{idx_dataset}, ...
    %         'MarkerFaceAlpha', 0.6, ...
    %         'MarkerEdgeAlpha', 0, ...
    %         'DisplayName', dataset_names{idx_dataset});
    % end
    %plot(snr_subject, acc_subject, '.m', 'DisplayName', 'Abs error', 'MarkerSize', 10)
    [correlation_error, pval] = corr(snr_current, bias_current);

    % p = polyfit(snr_current, bias_current, 1);
    % y_data = polyval(p, snr_current);
    % plot(snr_current, y_data, 'k', 'LineWidth', 2, 'DisplayName', 'Correlation')
    % if idx_approach == 1; ylabel('Bias [%]'); end
    % xlabel('SNR [db]')
    % ylim([-1 1])
    % %ylim([.45 .95]);
    % 
    % title(" \rho = " + num2str(round(correlation_error, 2)))

    fprintf('Correlation SNR and bias for %s: rho = %.3f, p = %.3f \n', ...
        get_fig_label_approach(approach_current), ...
        correlation_error, pval)
end

% l = legend('Interpreter', 'none');
% l.Layout.Tile = 'North';
% l.NumColumns = 3;

fontsize(gcf,scale=3)
exportgraphics(fig, 'REFTEP_preprocessing\Phase estimation scripts\figures\SNR.tif', 'Resolution', 500);

clear correlation_error correlation_bias p x y y_data y_lower y_upper ...
    idx_dataset idx_approach idx_subject data ans b subject_current


%% Methods figure

%for ind = [1 2 10 29 66 69]
ind = [69];
lim = [-300 100];
fig = figure(Units="normalized", Position=[0.5 0.4 0.25 0.5]);
t = tiledlayout(5, 2, 'TileSpacing','compact');
alpha = 1;

col = colors.blue;
col2 = colors.magenta;

% Ground truth
nexttile; hold on
title('Ground truth raw data', "FontName", 'TimesNewRoman')
ylabel('Ground truth', "FontName", 'TimesNewRoman')
plot(EEG_epoch_all.t_epoched_raw, EEG_epoch_all.data_epoched_raw(:,ind), 'Color', 'k', 'LineWidth', 2)
xlim(lim); xline(0, 'LineWidth',2)
xticks(0); yticks([]);
box on

% Filtered data
nexttile; hold on;
title('Ground truth filtered data', "FontName", 'TimesNewRoman')
data = EEG_epoch_all.data_cont_filt{1};
time = EEG_epoch_all.t_epoched_raw;
indices = round(EEG_epoch_all.events_cont{1}(69))+time-1;
plot(time, data(indices), 'Color', 'k', 'LineWidth', 2)
xlim(lim); xline(0, 'LineWidth',2); 
xticks(0); yticks([]);
yl = ylim;
ylim(yl)

ind_0 = round(EEG_epoch_all.t_data_filtered_PEAP) == 0;
box on

% PEAP: Unfiltered data with padding in other color
% Filtered data
% hilbert transform

% Raw data with interpolation
nexttile; hold on
title('Raw data, AR padding', "FontName", 'TimesNewRoman')
ylabel('PEAP', "FontName", 'TimesNewRoman')
ind_above_0 = EEG_epoch_all.t_data_filtered_PEAP >= -1;
plot(EEG_epoch_all.t_data_filtered_PEAP(ind_above_0), EEG_epoch_all.data_extended_PEAP(ind_above_0,ind), 'Color', [col alpha], 'LineWidth', 2)
plot(EEG_epoch_all.t_epoched_SSPE, EEG_epoch_all.data_epoched_SSPE(:,ind), 'Color', [0 0 0 alpha], 'LineWidth', 2)
xlim(lim); xline(0, 'LineWidth',2)
xticks(0); yticks([]);
box on

% Filtered data
nexttile; hold on;
title('Filtered data, AR padding', "FontName", 'TimesNewRoman')
plot(EEG_epoch_all.t_data_filtered_PEAP(~ind_above_0), EEG_epoch_all.data_filtered_PEAP(~ind_above_0,ind), 'Color', [0 0 0 alpha], 'LineWidth', 2)
plot(EEG_epoch_all.t_data_filtered_PEAP(ind_above_0), EEG_epoch_all.data_filtered_PEAP(ind_above_0,ind), 'Color', [col alpha], 'LineWidth', 2)
xlim(lim); xline(0, 'LineWidth',2); 
xticks(0); yticks([]);
ind_0 = round(EEG_epoch_all.t_data_filtered_PEAP) == 0;
yl = ylim;
ylim(yl)
box on

% PhastPadding
% Raw data with interpolation
nexttile; hold on
ylabel('PhastPadding', "FontName", 'TimesNewRoman')

title('Raw data, Phastimate padding', "FontName", 'TimesNewRoman')
ind_above_0 = EEG_epoch_all.t_data_filtered_PAR >= -1;
plot(EEG_epoch_all.t_data_filtered_PAR(ind_above_0), EEG_epoch_all.data_extended_PAR(ind_above_0,ind), 'Color', [col alpha], 'LineWidth', 2)
plot(EEG_epoch_all.t_epoched_SSPE, EEG_epoch_all.data_epoched_SSPE(:,ind), 'Color', [0 0 0 alpha], 'LineWidth', 2)
xlim(lim); xline(0, 'LineWidth',2)
xticks(0); yticks([]);
box on

% Filtered data
nexttile; hold on;
title('Filtered data, Phastimate padding', "FontName", 'TimesNewRoman')
plot(EEG_epoch_all.t_data_filtered_PAR(~ind_above_0), EEG_epoch_all.data_filtered_PAR(~ind_above_0,ind), 'Color', [0 0 0 alpha], 'LineWidth', 2)
plot(EEG_epoch_all.t_data_filtered_PAR(ind_above_0), EEG_epoch_all.data_filtered_PAR(ind_above_0,ind), 'Color', [col alpha], 'LineWidth', 2)
xlim(lim); xline(0, 'LineWidth',2)
ind_0 = round(EEG_epoch_all.t_data_filtered_PAR) == 0;
xticks(0); yticks([]);
yl = ylim;
ylim(yl)
box on

% Phastimate
% Raw data 
nexttile; hold on
ylabel('Phastimate', "FontName", 'TimesNewRoman')

title('Raw data, no padding', "FontName", 'TimesNewRoman')
plot(EEG_epoch_all.t_epoched_SSPE, EEG_epoch_all.data_epoched_SSPE(:,ind), 'Color', [0 0 0 alpha], 'LineWidth', 2)
xlim(lim); xline(0, 'LineWidth',2)
xticks(0); yticks([]);
box on

% Filtered data
nexttile; hold on;
title('Filtered data, Phastimate interpolation', "FontName", 'TimesNewRoman')

plot(EEG_epoch_all.t_data_filtered_default, EEG_epoch_all.data_filtered_default(:,ind), 'Color', [0 0 0 alpha], 'LineWidth', 2)
%plot(EEG_epoch_all.t_data_filtered_PAR(ind_above_0), EEG_epoch_all.data_filtered_PAR(ind_above_0,ind), 'Color', [col alpha], 'LineWidth', 2)
xlim(lim); xline(0, 'LineWidth',2)
% Phastimate interpolation
[~,~,~, data_temp] = phastimate_adapted(EEG_epoch_all.data_filtered_default(:,ind), 65, 30, 128);
time_temp = EEG_epoch_all.t_data_filtered_default(2:end) + 65;
ind_plot = time_temp > -65;
plot(time_temp(ind_plot), data_temp(ind_plot,:), 'Color', [col alpha], LineWidth=2)
yl = ylim;
ylim(yl)
yl = ylim;
ylim(yl)
rectangle('Position', [-65,  yl(1), 65, yl(2)-yl(1)], ...
          'EdgeColor', col2, 'Facecolor', col2, 'LineWidth', 2, 'FaceAlpha', 0.1)
ind_0 = round(time_temp) == 0;
xticks(0); yticks([]);
box on

% ETP
% Raw data 
nexttile; hold on
title('Raw data, no padding', "FontName", 'TimesNewRoman')
ylabel('ETP', "FontName", 'TimesNewRoman')

plot(EEG_epoch_all.t_epoched_SSPE, EEG_epoch_all.data_epoched_SSPE(:,ind), 'Color', [0 0 0 alpha], 'LineWidth', 2)
xlim(lim); xline(0, 'LineWidth',2)
xticks(0); yticks([]);
box on

% Interpolation
nexttile; hold on
xlim(lim); xline(0, 'LineWidth',2)
time = EEG_epoch_all.t_data_filtered_default;
data = EEG_epoch_all.data_filtered_default(:,ind);

% Find positive peaks
[pks, locs] = findpeaks(data, time);

% Keep only peaks before cutoff t = -40 ms
cutoff = -40;
validIdx = locs < cutoff;

% Last peak before cutoff
lastPeakTime = locs(find(validIdx, 1, 'last'));
lastPeakAmp  = pks(find(validIdx, 1, 'last'));
freq = 10;                  % Hz
omega = 2*pi*freq;
tInterp = lastPeakTime : 150;
t0 = lastPeakTime / 1000;    % seconds
tInterp_sec = tInterp / 1000;
interpSignal = lastPeakAmp * sin(omega * (tInterp_sec - t0) + pi/2);

plot(time, data, 'k', 'LineWidth', 2)
plot(tInterp, interpSignal, 'color', col, 'LineWidth', 2)
plot(lastPeakTime, lastPeakAmp, 'o', 'MarkerFaceColor', col)
yl = ylim;
ylim(yl)
rectangle('Position', [-40, yl(1), 40, yl(2)-yl(1)], ...
          'EdgeColor', col2, 'Facecolor', col2, 'LineWidth', 2, 'FaceAlpha', 0.1)
title('Filtered data, ETP interpolation', "FontName", 'TimesNewRoman')
ind_0 = round(tInterp) == 0;
title(t, 'Phase estimation methods', "FontName", 'TimesNewRoman')
xticks(0); yticks([]);
box on

%end
fontsize(gcf,scale=1.4)

% Export figure as tiff
exportgraphics(fig, 'REFTEP_preprocessing\Phase estimation scripts\figures\method.tif', 'Resolution', 500);


%% Calculate SNR of filtered data

SNR(EEG_epoch_all.data_cont_filt)

%% Calculate mean PSD

%% Quality control

fig = figure(Units="normalized", Position=[0.5 0.4 0.25 0.25]);
t = tiledlayout('flow', 'TileSpacing','compact');
title(t, 'Power-frequency spectrum', 'FontName', 'timesnewroman')

% testing set
% compute PSD for all subjects
clear snr_mu psd_all

for idx_participant = 1:length(EEG_epoch_all.t_cont)
    data = double(EEG_epoch_all.data_cont{idx_participant});
    [psd_temp, f] = pwelch(data,hamming(2*fs),[],[],fs);
    snr_mu.test(idx_participant) = 10*log10( ...
        bandpower(psd_temp,f,[9 13],'psd') / ...
        (bandpower(psd_temp,f,[0 fs/2],'psd') - bandpower(psd_temp,f,[9 13],'psd')) );
    psd_all.test(:,idx_participant) = psd_temp(f<100);
    data = double(EEG_epoch_all.data_cont_filt{idx_participant});
    [psd_temp, f] = pwelch(data,hamming(2*fs),[],[],fs);
    snr_mu.test_filt(idx_participant) = 10*log10( ...
        bandpower(psd_temp,f,[9 13],'psd') / ...
        (bandpower(psd_temp,f,[0 fs/2],'psd') - bandpower(psd_temp,f,[9 13],'psd')) );
end

outliercount.test = sum(EEG_epoch_all.n_rejection);

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

% outlier count
fprintf('Rejected epochs (total): %d\n (relative): %.2f percent\n', outliercount.test, outliercount.test/length(EEG_epoch_all.subject)*100);

% Plot PSD
f_plot   = f(f > 1 & f < 45);
psd_plot = log(psd_all.test((f > 1 & f < 45), :));

nexttile; hold on; grid on

% plot(f_plot, mean(psd_plot,2), 'Color', 'w', 'LineWidth', 2, 'DisplayName', '\bf Test')
% hold on

plot(f_plot, mean(psd_plot,2), 'Color', colors.cyan, 'LineWidth', 2, 'DisplayName', 'Test')

% Plot per subgroup
% dataset_plot_names = ["Young control", "Elderly control" "Stroke intact" "Stroke affected"];
% linestyle_list = ["-", "--", "-.", ":"];
% hold on
% idx_dataset = floor(subject_dataset_labels/100);
% idx_dataset(idx_dataset == 1) = 2;
% idx_dataset = idx_dataset - 1;
% for i = 1:4
%     psd_subplot = psd_plot(:,idx_dataset == i);
%     plot(f_plot, mean(psd_subplot,2), 'Color', colors.cyan, 'LineWidth', 1, 'DisplayName', dataset_plot_names{i}, 'Linestyle', linestyle_list(i))
% end

% train set
filepath1 = [path.load 'ground_truth_' 'train'];
load([filepath1 '\EEG_phase_ground_truth_' 'train' '_total']);

for idx_participant = 1:length(EEG_epoch_all.t_cont)
    data = double(EEG_epoch_all.data_cont{idx_participant});
    [psd_temp, f] = pwelch(data,hamming(2*fs),[],[],fs);
    snr_mu.train(idx_participant) = 10*log10( ...
        bandpower(psd_temp,f,[9 13],'psd') / ...
        (bandpower(psd_temp,f,[0 fs/2],'psd') - bandpower(psd_temp,f,[9 13],'psd')) );
    psd_all.train(:,idx_participant) = psd_temp(f<100);
    data = double(EEG_epoch_all.data_cont_filt{idx_participant});
    [psd_temp, f] = pwelch(data,hamming(2*fs),[],[],fs);
    snr_mu.train_filt(idx_participant) = 10*log10( ...
        bandpower(psd_temp,f,[9 13],'psd') / ...
        (bandpower(psd_temp,f,[0 fs/2],'psd') - bandpower(psd_temp,f,[9 13],'psd')) );
end

outliercount.train = sum(EEG_epoch_all.n_rejection);

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

% outlier count
fprintf('Rejected epochs (total): %d\n (relative): %.2f percent\n', outliercount.train, outliercount.train/length(EEG_epoch_all.subject)*100);

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
load([filepath '\phase_' type '_total'])
fontsize(gcf,scale=1.4)
set(gca,'FontName', 'Times')

% Plot per subgroup
% figure
% linestyle_list = ["-", "--", "-.", ":"];
% hold on
% idx_dataset = floor(subject_dataset_labels/100);
% idx_dataset(idx_dataset == 1) = 2;
% idx_dataset = idx_dataset - 1;
% for i = 1:4
%     psd_subplot = psd_plot(:,idx_dataset == i);
%     plot(f_plot, mean(psd_subplot,2), 'Color', colors.magenta, 'LineWidth', 1, 'DisplayName', dataset_plot_names{i}, 'Linestyle', linestyle_list(i))
% end

exportgraphics(fig, 'REFTEP_preprocessing\Phase estimation scripts\figures\PSD_all.tif', 'Resolution', 500);



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

