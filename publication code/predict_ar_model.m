function [data_extended, time_extended, time_measure] = predict_ar_model(data, time, order_AR, n_output, approach_AR)
%PREDICT_AR Performs AR-based forward prediction using the specified approach
%
% Inputs:
%   - data:        [n_samples x n_epochs] matrix of input data
%   - order_AR:    scalar AR model order
%   - n_output:    scalar number of values to predict
%   - approach_AR: string, AR model estimation method
%
% Outputs:
%   - data_extended:       [n_samples + n_output_current x n_epochs] matrix with predictions
%   - time_measure:        scalar, computation time in seconds

% Prepare output matrix
data_extended = [data; zeros(n_output, size(data, 2))];
time_extended = [time, 0:1:n_output-1];


switch approach_AR
    case "burg"
        tic
        model = arburg(data, order_AR);
        coefficients = -1 * flip(model(:, 2:end)');
        for i = n_output:-1:1
            data_extended(end-i+1,:) = ...
                sum(coefficients .* data_extended((end-i-order_AR+1):(end-i),:));
        end
        time_measure = toc;

    case "yw"
        tic
        model = aryule(data, order_AR);
        coefficients = -1 * flip(model(:, 2:end)');
        for i = n_output:-1:1
            data_extended(end-i+1,:) = ...
                sum(coefficients .* data_extended((end-i-order_AR+1):(end-i),:));
        end
        time_measure = toc;

    case "arcov"
        tic
        model = arcov(data, order_AR);
        coefficients = -1 * flip(model(:, 2:end)');
        for i = n_output:-1:1
            data_extended(end-i+1,:) = ...
                sum(coefficients .* data_extended((end-i-order_AR+1):(end-i),:));
        end
        time_measure = toc;

    case "armcov"
        tic
        model = armcov(data, order_AR);
        coefficients = -1 * flip(model(:, 2:end)');
        for i = n_output:-1:1
            data_extended(end-i+1,:) = ...
                sum(coefficients .* data_extended((end-i-order_AR+1):(end-i),:));
        end
        time_measure = toc;

    case "lpc"
        tic
        model = lpc(double(data), order_AR);
        coefficients = -1 * flip(model(:,2:end))';  % Make it a row vector
        for i = n_output:-1:1
            data_extended(end - i + 1, :) = ...
                sum(coefficients .* data_extended((end - i - order_AR + 1):(end - i), :));
        end
        time_measure = toc;

    otherwise
        tic
        n_epochs = size(data, 2);
        n_samples = size(data, 1);
        for idx_sample = 1:n_epochs
            if mod(idx_sample, 1000) == 0
                fprintf('sample %d of %d\n', idx_sample, n_epochs);
            end
            sample = double(data(:, idx_sample));
            model = ar(sample, order_AR, 'Approach', approach_AR);
            past_values = sample(end-order_AR+1:end);
            predicted = zeros(n_output, 1);
            for i = 1:n_output
                next_val = -model.A(2:end) * flipud(past_values);
                predicted(i) = next_val;
                past_values = [past_values(2:end); next_val];
            end
            data_extended(n_samples+1:end, idx_sample) = predicted;
        end
        time_measure = toc;
end
end
