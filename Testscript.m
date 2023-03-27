% Test script for SpatioTemporalInnovationsFilterDesigner, which estimates an innovations filter from ensembles of multichannel random
% processes.
%
% The Open-Source Electrophysiological Toolbox (OSET)
% https://github.com/alphanumericslab/OSET
% Reza Sameni, Feb 2023

close all
clear
clc

case_study = 'EEG'; % 'EEG' or 'PCG'
switch case_study

    case 'EEG' % A sample EEG from OSET
        load 'Desktop/Patient_1_interictal_segment_0001.mat';
        fs = 5000; % Sampling frequency (5000)
        num_segments = 1; % Number of segments used to split the large data matrix (trying with one first)
        wlen = round(fs * 600.0); % length of each segment (600 sec)
        sample_data = cell(1, num_segments);
        for k = 1 : num_segments
            sample_data{k} = interictal_segment_1.data(:, (k-1)*wlen + 1 : k*wlen);
        end
    otherwise
        error('Undefined case study. Try EEG, PCG, or adapt the script according to your data.')
end

params.spatial_filter_type = 'PCA'; % 'BY_PASS', 'PCA' or 'ICA'
params.normalize_records = true; % true/false
params.fs = fs; % sampling frequency
params.keep_mean = true; % true/false
params.spectral_len = 512; % number of frequency bins for spectral estimation (513?_
params.filter_len = params.spectral_len; % innovations filter length (best practice to set equal to params.spectral_len)
params.smooth_spectrum = false; % true/false smooth the spectrum before spectral factorization
if params.smooth_spectrum
    params.lambda = 10000.0; % Tikhonov regularization factor used for spectral smoothing
end
params.spectral_averaging_method = 'MEDIAN'; % 'MEDIAN', 'MEAN', 'MAX', 'MIN', 'MAX_MIN_AVG'
params.innovation_filter_type = 'LINEAR_PHASE'; % 'LINEAR_PHASE', 'MIN_PHASE'; (linear for EEG waves)
params.plot_results = true; % true/false plot results

% Design the innovations filter
[h_innovations, A, S_mean, S_median, S_max, S_min, S_max_min_avg] = SpatioTemporalInnovationsFilterDesigner(sample_data, params);


% Test the innovations filter
dat = sample_data{1};
syn_signal_len = size(dat, 2);
x = randn(length(h_innovations), syn_signal_len); % generating white noise
y = randn(length(h_innovations), syn_signal_len); % synth data stored in y
for ch = 1 : length(h_innovations) % loop iterates through each channel of the EEG signal
    y(ch, :) = filter(h_innovations{ch}, 1, x(ch, :));
    [H,F] = freqz(h_innovations{ch}, 1, [], fs);

    figure
    subplot(411)
    plot((0 : length(dat(ch, :)) - 1)/fs, dat(ch, :));
    grid
    xlabel('time(s)');
    ylabel('real data, amplitude');
    set(gca, 'fontsize', 18)

    subplot(412)
    plot((0 : length(y(ch, :)) - 1)/fs, y(ch, :));
    grid
    xlabel('time(s)');
    ylabel('innovations process, amplitude');
    set(gca, 'fontsize', 18)

    subplot(413)
    plot(F, 20*log10(abs(H)));
    grid
    xlabel('frequency(Hz)');
    ylabel('innovations filter magnitude response');
    set(gca, 'fontsize', 18)

    subplot(414)
    plot(F, unwrap(angle(H)));
    grid
    xlabel('frequency(Hz)');
    ylabel('innovations filter phase response');
    set(gca, 'fontsize', 18)
    sgtitle('Real data vs synthetic innovations process', 'fontsize', 18);
end

% save InnovationsFilterDesign