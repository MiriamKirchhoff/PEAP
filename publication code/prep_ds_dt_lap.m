function EEG_lap = prep_ds_dt_lap(EEG, hjorth)

if nargin < 2
    hjorth.channel = 'C3';
    hjorth.electrodes = {'FCC5h', 'FCC3h', 'CCP5h', 'CCP3h'};
end


%% Settings

temp = struct2cell(EEG.chanlocs);
hjorth.channel_idx = find(strcmpi(hjorth.channel, temp(1,:)));
for i = 1:length(hjorth.electrodes)
    hjorth.electrode_idx(i) = find(strcmpi(hjorth.electrodes{i}, temp(1,:)));
end
clear i


%% Downsample

% settings
ds_fs = 1000;                               % downsampling freqency in Hz

% start timer
tic

% start notification
disp('Start downsampling...')

% downsample data
EEG = pop_resample(EEG, ds_fs);


%% Detrend

% settings
% define prestim period
t.detrend = [min(EEG.times), max(EEG.times)];         % ms time wrt TMS

% start timer
tic

% start notification
disp('Start detrending...')

% detrend data linearly
EEG = tesa_detrend(EEG, 'linear', t.detrend);

% timing notification
fprintf('detrending completed. Expired time: %.0f seconds \n', toc)


%% Laplacian montage

% start timer
tic

% start notification
disp('Start laplacian...')

EEG_lap = EEG;

% Sanity check: EEG.data should be [nbchan x pnts x trials] or [nbchan x pnts]
data = EEG_lap.data;

% If 2D, reshape to 3D for consistent processing
if ndims(data) == 2
    data = reshape(data, size(data,1), size(data,2), 1);
end

% apply laplacian montage
lap_data = data(hjorth.channel_idx,:,:) - ...
    mean(data(hjorth.electrode_idx,:,:), 1);

clear data

% Update EEG_lap structure
EEG_lap.data = lap_data;
EEG_lap.nbchan = 1;
EEG_lap.chanlocs = [hjorth.channel '-centered laplacian'];
EEG_lap.soundLeadfield = [];

% timing notification
fprintf('Laplacian montage completed. Expired time: %.0f seconds \n', toc)


end % eof [EEG_lap] = preprocess(EEG)