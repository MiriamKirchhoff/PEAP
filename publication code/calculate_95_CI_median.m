function [med, CI_lower, CI_upper] = calculate_95_CI_median(data, nBoot)
% calculate_95_CI_median computes the row-wise median and 95% CI using bootstrapping
%
% INPUT:
%   data  - an nRows x nSamples matrix
%   nBoot - number of bootstrap resamples (default = 1000)
%
% OUTPUTS:
%   med       - median across rows
%   CI_lower  - lower bound of 95% confidence interval
%   CI_upper  - upper bound of 95% confidence interval

if nargin < 2
    nBoot = 1000;
end

nCols = size(data, 2);

med = median(data, 1, 'omitnan');
CI_lower = nan(1, nCols);
CI_upper = nan(1, nCols);

for j = 1:nCols
    x = data(:, j);
    x = x(~isnan(x));
    if numel(x) < 2
        med(j) = NaN;
        CI_lower(j) = NaN;
        CI_upper(j) = NaN;
        continue
    end

    bootMedians = bootstrp(nBoot, @median, x);
    CI = prctile(bootMedians, [2.5 97.5]);
    CI_lower(j) = CI(1);
    CI_upper(j) = CI(2);
end
end