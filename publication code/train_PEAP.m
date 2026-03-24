function ar_parameters_optimal = train_PEAP
% function ar_parameters_optimal = train_PEAP
% Function to get optimal PEAP parameters.
%
% OUTPUT 
% ar_parameters_optimal: struct with fields
    % n_input: Number of input samples
    % n_output: Number of output samples
    % order: AR model order
    % approach: AR model used (burg or yule-walter (yw))


%% Settings

path.load = '';
path.save = '';

% Always use training data
type = 'train';


%% Load ground truth and epoched data

filepath = [path.save 'ground_truth_' type];
load([filepath '\EEG_phase_ground_truth_' type '_total']);
load('settings.mat');

%% Tune parameters: n_input, n_output, AR order, AR approach

data = EEG_epoch_all.data_epoched_SSPE;

n_input = 700:20:1000;
n_output = 110:20:300;  % round(linspace(100, 150, 3));
order_AR = 110:10:200; % round(linspace(75, 150, 4));
approach_AR =  ["burg", "yw"]; %["burg", "yw", "fb", "gl", "ls", "arcov", "armcov"];
optimization_time = -1;      % [ms]

n_combinations = length(n_input) * length(n_output) * length(order_AR) * length(approach_AR);


combination_idx  = 1;
clear rmse bias acc_relevant
rmse = [];
bias = [];
acc_relevant = [];
% figure; hold on; grid on;

% kuiper sensitivity: indicates min pval to be reached by kupier test.
% pval indicates
kuiper_threshold = .2;

for idx_method = 1:length(approach_AR)
    for idx_order = 1:length(order_AR)
        for idx_in = 1:length(n_input)

            % Get current parameters
            n_input_current = n_input(idx_in);
            order_AR_current = order_AR(idx_order);
            approach_AR_current = approach_AR(idx_method);
            % get max output to already calculate all points and sub-select
            % later
            n_output_current = max(n_output);
            
            % get as many samples prior to 0 as indicated
            data_current = data(end-n_input_current+1:end, :);
            time_current = EEG_epoch_all.t_epoched_SSPE(end-n_input_current+1:end);


            % Check that AR order < input length to avoid errors
            if order_AR_current >= n_input_current
                fprintf('Skipping: AR order (%d) >= input length (%d)\n', order_AR_current, n_input_current);
                rmse = NaN;
                continue;
            end

            %% Forecast AR using specified method
            [data_extended, time_extended] = predict_ar_model( ...
                data_current, time_current, order_AR_current, ...
                n_output_current, approach_AR_current);


            for idx_out = 1:length(n_output)

                combination_idx = ...
                    (idx_method - 1) * (length(order_AR) * length(n_input) * length(n_output)) + ...
                    (idx_order  - 1) * (length(n_input) * length(n_output)) + ...
                    (idx_in     - 1) * (length(n_output)) + ...
                    idx_out;

                % Get current parameters
                n_output_current = n_output(idx_out);


                fprintf('Combination %4d of %4d: method=%s | order=%d | input=%d | output=%d | \n', ...
                    combination_idx, n_combinations, approach_AR_current, order_AR_current, n_input_current, n_output_current);

                % subselect data until output size is reached
                data_extended_current = data_extended(1:n_input_current+n_output_current,:);
                time_extended_current = time_extended(1:n_input_current+n_output_current);


                %% Filter data

                data_extended_current = filtfilt(settings.filter, 1, data_extended_current);


                %% Calculate phase

                hilb = hilbert(data_extended_current);

                phase_hilb = angle(hilb(round(time_extended_current)==optimization_time,:));


                %% Calculate accuracy

                reference_phase_location = find(round(phase_full_all.time) == optimization_time);
                gt_current = phase_full_all.ground_truth(:, reference_phase_location)';

                rmse(combination_idx) = sqrt(mean(circ_dist(phase_hilb, gt_current).^2));
                bias(combination_idx) = mean(circ_dist(phase_hilb, gt_current));

                % make a matrix with all locations for the full gt
                % TODO
                % locs_epochs = repmat(epochs_current', 1, numel(win_forecast));
                % locs_epochs = locs_epochs + repmat(win_forecast, numel(epochs_current), 1);
                % phase_full.ground_truth = angle(hilb(locs_epochs));


                %% Save phase? Or just accuaracy?

                disp(['rmse = ' num2str(rmse(combination_idx))])

                if isempty(acc_relevant)
                    
                    pval = circ_kuipertest(phase_hilb, gt_current, 1000, 0);
                    disp(['kupier pval = ' num2str(pval)])
                    if pval > kuiper_threshold
                        acc_relevant(combination_idx) = rmse(combination_idx);
                        disp('kuiper test passed')
                        disp(['rmse = ' num2str(rmse(combination_idx))])

                        if min(acc_relevant) == acc_relevant(combination_idx)

                            title_current = strjoin(['Combination ' num2str(combination_idx) ...
                                ': method=' approach_AR_current ...
                                ' | order=' num2str(order_AR_current) ...
                                ' | input=' num2str(n_input_current) ...
                                ' | output=' num2str(n_output_current) ...
                                ' | rmse=' num2str(rmse(combination_idx))]);
                            %%
                            plot_stats(phase_hilb, rmse, bias, gt_current, EEG_epoch_all, title_current)
                            %%
                            drawnow

                            ar_parameters_optimal.n_input = n_input_current;
                            ar_parameters_optimal.n_output = n_output_current;
                            ar_parameters_optimal.order = order_AR_current;
                            ar_parameters_optimal.approach = approach_AR_current;

                        end
                    else
                        acc_relevant(combination_idx) = NaN;
                    end

                    
                % if there is any improvement
                % elseif rmse(combination_idx) < min(min(acc_relevant), pi)

                % if there is an improvement of at least 0.1% w.r.t. prior
                % accuracy
                elseif min(min(acc_relevant), pi) - rmse(combination_idx) > min(min(acc_relevant), pi)*0.001
                    pval = circ_kuipertest(phase_hilb, gt_current, 1000, 0);
                    disp(['kupier pval = ' num2str(pval)])
                    if pval > kuiper_threshold
                        acc_relevant(combination_idx) = rmse(combination_idx);
                        disp('kuiper test passed')
                        

                        if min(acc_relevant) == acc_relevant(combination_idx)

                            title_current = strjoin(['Combination ' num2str(combination_idx) ...
                                ': method=' approach_AR_current ...
                                ' | order=' num2str(order_AR_current) ...
                                ' | input=' num2str(n_input_current) ...
                                ' | output=' num2str(n_output_current) ...
                                ' | rmse=' num2str(rmse(combination_idx))]);
                            plot_stats(phase_hilb, rmse, bias, gt_current, EEG_epoch_all, title_current)
                            drawnow

                            ar_parameters_optimal.n_input = n_input_current;
                            ar_parameters_optimal.n_output = n_output_current;
                            ar_parameters_optimal.order = order_AR_current;
                            ar_parameters_optimal.approach = approach_AR_current;

                        end
                    else
                        acc_relevant(combination_idx) = NaN;
                    end

                else
                    acc_relevant(combination_idx) = NaN;
                end

            end
        end
    end
    % plot at end of each predictor
    title_current = strjoin(['Combination ' num2str(combination_idx) ...
        ': method=' approach_AR_current ...
        ' | order=' num2str(order_AR_current) ...
        ' | input=' num2str(n_input_current) ...
        ' | output=' num2str(n_output_current) ...
        ' | rmse=' num2str(rmse(combination_idx))]);
    plot_stats(phase_hilb, rmse, bias, gt_current, EEG_epoch_all, title_current)
    drawnow
end



%% Save data

filepath = [path.save '\algorithm_training_results\'];

if ~exist(filepath, 'dir')
    mkdir(filepath)
end

save([filepath 'PEAP'], 'ar_parameters_optimal', '-v7.3')






