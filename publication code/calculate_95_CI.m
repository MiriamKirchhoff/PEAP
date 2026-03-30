function [mu, CI_lower, CI_upper] = calculate_95_CI(data)
% rowMeanCI95 computes the row-wise mean and 95% confidence interval
%
% INPUT:
%   data - an nRows x nSamples matrix
%
% OUTPUTS:
%   mu       - mean across columns (per row)
%   CI_lower - lower bound of 95% confidence interval
%   CI_upper - upper bound of 95% confidence interval

% Number of samples (columns)
n = size(data, 1);

% Row-wise mean
mu = mean(data, 'omitnan');

% Row-wise standard deviation (normalized by n-1)
sigma = std(data, 0, 'omitnan');

% Standard error of the mean
SEM = sigma / sqrt(n);

% t-value for 95% CI, two-tailed
t_val = tinv(0.975, n - 1);

% Confidence intervals
CI_lower = mu - t_val * SEM;
CI_upper = mu + t_val * SEM;
end