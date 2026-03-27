function phastimate_parameters_optimal = train_phastimate
% this is a subject-specific version of Phastimate individualized setting 
% search (Zrenner 2018, 2020)

%% Train phastimate settings based on resting state

path.load = '';
path.save = '';

type = 'train';


%% Load ground truth and epoched data

filepath = [path.save 'ground_truth_' type];
load([filepath '\EEG_phase_ground_truth_' type '_total'])


%% Settings

subjects = unique(EEG_epoch_all.subject);


for idx_subject = 1:length(subjects)  % iterate over subjects

    subject_current = subjects(idx_subject);
    disp(['start ' num2str(idx_subject) ' of ' num2str(length(subjects))])
    data = EEG_epoch_all.data_cont{idx_subject};
    fs = EEG_epoch_all.fs;


    %% Estimate phase using phastimate

    %% preliminaries

    % check for toolboxes
    assert(~isempty(which('designfilt')), 'filter design function designfilt.m not found, is the Signal Processing Toolbox installed?')
    assert(~isempty(which('range')), 'statistical function range.m not found, is the Statistics and Machine Learning Toolbox installed?')
    assert(~isempty(which('ga')), 'genetic algorithm function ga.m not found, is the Global Optimization Toolbox installed?')

    % circular statistics functions (simplified from circstat toolbox)
    ang_mean = @(x) angle(mean(exp(1i*x)));
    ang_diff = @(x, y) angle(exp(1i*x)./exp(1i*y));
    ang_var = @(x) 1-abs(mean(exp(1i*x)));
    %ang_var2dev = @(v) sqrt(2*v); % circstat preferred formula uses angular deviation (bounded from 0 to sqrt(2)) which is sqrt(2*(1-r))
    ang_var2dev = @(v) sqrt(-2*log(1-v)); % formula for circular standard deviation is sqrt(-2*ln(r))

    %% set constants

    % filter design method for phastimate (order and peak frequency is variable)
    design_phastimate_filter = @(ord, freq, fs) designfilt('bandpassfir', 'FilterOrder', ord, 'CutoffFrequency1', freq-1, 'CutoffFrequency2', freq+1, 'SampleRate', fs, 'DesignMethod', 'window');

    NUM_EPOCHS = 497;
    HILBERTWIN = 128; % this is an appropriate window for alpha at 1000 Hz
    PEAK_FREQUENCY_INTERVAL = [8 14];

    %% load resting sate data into a master data table 'T'

    % only select first three minutes
    time = 1:fs*60*3;
    T = table();
    T.data = data;
    T.fs = fs * ones(height(T),1);


    %% Step 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine alpha peak frequency and signal to noise ratio
    % Note:
    % - data is not cleaned before this step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    T.peak_frequency = nan(height(T),1);
    T.peak_SNR = nan(height(T),1);

    epochs = create_epochs_overlapping(T.data, T.fs); % from continuous data

    T.peak_frequency = [];

    % estimate SNR and plot
    [peak_frequency, peak_SNR] = estimate_SNR(epochs, T.fs, PEAK_FREQUENCY_INTERVAL, []);
    if isempty(peak_frequency)
        peak_frequency = 10;
            peak_SNR = NaN;
    end

    T.peak_frequency = peak_frequency;
    T.peak_SNR = peak_SNR;
    
    clear('subplot_index', 'row_index',  'ax', 'epochs', 'peak_frequency', 'peak_SNR', 'i')


    %% Step 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % determine true phase and amplitude, as well as variance of "true" phase
    % - epochs are recreated from the data when needed to save memory space
    %
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % add columns for epoch-by-epoch data
    T.epochs_truephase_mean = nan(height(T), NUM_EPOCHS);
    T.epochs_truephase_ang_var = nan(height(T), NUM_EPOCHS);
    T.epochs_trueamp_mean = nan(height(T), NUM_EPOCHS);
    T.epochs_trueamp_var = nan(height(T), NUM_EPOCHS);

    epochs = create_epochs_overlapping(T.data, T.fs); % from continuous data

    peak_frequency = double(T.peak_frequency);

    % set-up equivalent filter objects for given peak frequency
    filter_objects = {};
    fs = T.fs;
    for ord = [2 3 4 5] % FIR - windowed sinc
        filter_objects = {filter_objects{:} designfilt('bandpassfir', 'FilterOrder', round(ord * (fs/peak_frequency)), 'CutoffFrequency1', peak_frequency-1, 'CutoffFrequency2', peak_frequency+1, 'SampleRate', fs, 'DesignMethod', 'window')};
    end
    for ord = [3 4 5] % FIR - least squares (equiripple is similar)
        filter_objects = {filter_objects{:} designfilt('bandpassfir', 'FilterOrder', round(ord * (fs/peak_frequency)), 'StopbandFrequency1', peak_frequency-4, 'PassbandFrequency1', peak_frequency-1, 'PassbandFrequency2', peak_frequency+1, 'StopbandFrequency2', peak_frequency+4, 'SampleRate', fs, 'DesignMethod', 'ls')};
    end
    for ord = [4 8 12] % IIR - butterworth
        filter_objects = {filter_objects{:} designfilt('bandpassiir', 'FilterOrder', ord, 'HalfPowerFrequency1', peak_frequency-1, 'HalfPowerFrequency2', peak_frequency+1, 'SampleRate', fs, 'DesignMethod', 'butter')};
    end
    for ord = [4 6 8] % IIR - chebychev I
        filter_objects = {filter_objects{:} designfilt('bandpassiir', 'FilterOrder', ord, 'PassbandFrequency1', peak_frequency-1.5, 'PassbandFrequency2', peak_frequency+1.5, 'PassbandRipple', 0.5, 'SampleRate', fs, 'DesignMethod', 'cheby1')};
    end
    for attenuation = [10 20] % IIR - elliptic
        filter_objects = {filter_objects{:} designfilt('bandpassiir', 'StopbandFrequency1', peak_frequency-2, 'PassbandFrequency1', peak_frequency-1, 'PassbandFrequency2', peak_frequency+1, 'StopbandFrequency2', peak_frequency+2, 'StopbandAttenuation1', attenuation, 'PassbandRipple', 0.5, 'StopbandAttenuation2', attenuation, 'SampleRate', fs, 'DesignMethod', 'ellip', 'MatchExactly', 'passband')};
    end

    [truephase_mean, truephase_variance, trueamp_mean, trueamp_variance] = phastimate_truephase(epochs, filter_objects);

    T.epochs_truephase_mean = truephase_mean;
    T.epochs_truephase_ang_var = truephase_variance;

    T.epochs_trueamp_mean = trueamp_mean;
    T.epochs_trueamp_var = trueamp_variance;

    clear('subplot_index', 'row_index', 'epochs', 'peak_frequency', 'filter_objects', 'fs', 'ord', 'truephase_mean', 'truephase_variance', 'trueamp_mean', 'trueamp_variance', 'ax', 'h', 'i')


    %% Step 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine optimized phastimate parameters and resulting estimate
    %
    %
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %add relevant columns
    T.optim_window = nan(height(T), 1);
    T.optim_filter_ord = nan(height(T), 1);
    T.optim_edge = nan(height(T), 1);
    T.optim_ar_ord = nan(height(T), 1);
    T.optim_fval = nan(height(T), 1);

    fprintf('\nNow running genetic algorithm to find optimized phastimate parameters...')

    epochs = create_epochs_overlapping(T.data, T.fs); % from continuous data

    peak_frequency = double(T.peak_frequency);

    filter_order_range = 100:250;

    filter_objects_by_order = {}; %the index has to correspond to the order for the genetic algorithm
    for ord = filter_order_range
        filter_objects_by_order{ord} = design_phastimate_filter(ord, peak_frequency, T.fs);
    end

    bounds_filter_order = [filter_order_range(1) filter_order_range(end)];
    bounds_window = [400 750];
    bounds_edge = [30 120];
    bounds_ar_order = [5 60];

    % the includemask allows optimizing for a subset of epochs
    % it makes sense to exclude epochs that would also be excluded by the
    % real-time system, e.g. if artifacts are detected so as to not optimize
    % for noisy epochs that wouldn't result in a stimulus anyway

    % subselect according to truephase variance
    %includemask = T.epochs_truephase_angdev <= quantile(T.epochs_truephase_angdev, 0.5);

    % subselect according to true amplitude
    includemask = T.epochs_trueamp_mean >= quantile(T.epochs_trueamp_mean, 0.5);

    [optimal_parameters, ga_output] = phastimate_optimize(epochs(1:ceil(end/2),includemask), T.epochs_truephase_mean(includemask), filter_objects_by_order, bounds_filter_order, bounds_window, bounds_edge, bounds_ar_order, HILBERTWIN);

    % rerun phastimate with the optimized settings to confirm result
    D = design_phastimate_filter(optimal_parameters.filter_order, peak_frequency, T.fs);
    [estphase, estamp] = phastimate(epochs(((-optimal_parameters.window_length+1):0)+ceil(end/2),:), D, optimal_parameters.edge, optimal_parameters.ar_order, 128);

    % sanity check if the angular deviation matches the result of the optimization
    phases_error = ang_diff(T.epochs_truephase_mean, estphase);
    assert(abs(optimal_parameters.fval - ang_var(phases_error(includemask))) < 0.01, 'could not replicate result of optimization, were the same filters used?')

    T.optim_window = optimal_parameters.window_length;
    T.optim_filter_ord = optimal_parameters.filter_order;
    T.optim_edge = optimal_parameters.edge;
    T.optim_ar_ord = optimal_parameters.ar_order;
    T.optim_fval = optimal_parameters.fval;

    fprintf('\nDone.\n')

    clear('row_index',  'epochs', 'peak_frequency', 'filter_order_range', 'ord', 'bounds_ar_order', 'filter_objects_by_order', 'bounds_edge', 'bounds_filter_order', 'bounds_window', 'D', 'includemask', 'estamp', 'estphase', 'optimal_parameters', 'ga_output', 'phases_error')

    %% Save data in array

    phastimate_parameters_optimal.settings{idx_subject} = T;
    phastimate_parameters_optimal.subject{idx_subject} = subject_current;


end % idx_subject = subjects  % iterate over subjects

%% Save data

filepath = [path.save '\algorithm_training_results\'];

if ~exist(filepath, 'dir')
    mkdir(filepath)
end

save([filepath 'phastimate'], 'phastimate_parameters_optimal', '-v7.3')


end
