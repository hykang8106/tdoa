function [position_error_torrieri, x_est_torrieri, cep, gdop, lambda1, lambda2, theta_rad] = ...
    sub_simulate_tdoa_sensor_fixed(sensor_position, snr_db, target_position, ...
    plot_signal, plot_position, use_only_torrieri_method, target_radius_meter, rician_param, ...
    tx_signal, fs, nfft, bw_mhz)
% subroutine to simulate target location using tdoa method
%
% [input]
% - sensor_position: sensor position, dimension = sensor_length x 2
% - snr_db: snr in db
% - target_position: target position, dimension = 1 x 2 
% - plot_signal: boolean. to control plot of tx, rx, correlation
% - plot_position: boolean. to control plot of target, sensor, tdoa curve
% - use_only_torrieri_method: boolean. 
%   0 = use all method(least square estimator, linear equation solver, hyperbolic equation solver)
%   1 = use least square estimator
% - target_radius_meter: uca_radius_meter * radius_ratio
% - rician_param: rician fading parameter
% - tx_signal: target signal(lte prs)
% - fs: target signal(lte prs) sample rate
% - nfft: subcarrier number in lte prs
% - bw_mhz: target signal(lte prs) bandwidth in mhz
%

sensor_length = size(sensor_position, 1);

%% Plot Location of rf sensor and target emitter 

if plot_position
    H = plot_sensor_position(sensor_position, target_position, target_radius_meter);
end

%% Transmitter(target emitter: lte base station) Configuration

tx_signal_length = length(tx_signal);

% Plot transmit waveforms from target emitter
if plot_signal
    plot_tx_signal(tx_signal, fs);
end
fprintf('sample rate = %f MHz, distance per sample = %f meter\n', ...
    fs / 1e6, physconst('LightSpeed') / fs);

%% Compute delays from target emitter to rf sensors

speedOfLight = physconst('LightSpeed'); % Speed of light in m/s

sampleDelay = zeros(1, sensor_length);
radius = zeros(1, sensor_length);
% radius = cell(1, sensor_length);
faded_signal = zeros(tx_signal_length, sensor_length);
for i = 1 : sensor_length
    [~, radius(i)] = cart2pol(sensor_position(i, 1) - target_position(1), ...
        sensor_position(i, 2) - target_position(2));
    delay = radius(i) / speedOfLight;                 % Delay in seconds
    
    if ~isempty(rician_param)
        k = rician_param{i, 1};
        if ~isempty(k)
            % tau in rician_param is delta tau ratio(not absolute delay), see rician param excel file
            tau = rician_param{i, 2} * delay;
            % ########## For outdoor environments,
            % ########## path delays after the first are typically between 100 ns(= 1e-7 s) and 10 us(= 1e-5 s)
            % ########## 100 ns = 30 m, 10 us = 3 km
            if sum(tau(2:end) > 1e-5) || sum(tau(2:end) < 1e-7)
                tau
                fprintf('### [warning] path delay must be 1e-7 ~ 1e-5 sec\n');
            end
            pdb = rician_param{i, 3};
            [faded_signal(:, i), chan] = rician_fading(tx_signal, fs, k, tau, pdb);
        else
            faded_signal(:, i) = tx_signal;
        end
    end
    
%     sampleDelay(i) = floor(delay * fs); % Delay in samples
%     sampleDelay(i) = ceil(delay * fs); % Delay in samples
%     sampleDelay(i) = fix(delay * fs); % Delay in samples
    sampleDelay(i) = round(delay * fs); % Delay in samples
end
sampleDelay

%% Create Received Waveforms in each rf sensor and Plot Received Waveform

fc = 2.1e9; % lte prs carrier freq

% ###########################################################################
% ## for carrier freq of fsk 422mhz signal, lte carrier freq will be used
% ## because path loss difference is very small:
% ## hPositioningPathLoss(3e3, 2.1e9) => pathloss = 128 db
% ## hPositioningPathLoss(3e3, 422e6) => pathloss = 127 db
% ###########################################################################

rx_signal = zeros(tx_signal_length + max(sampleDelay), sensor_length);
for i = 1 : sensor_length
    % Urban Macro LOS path loss per TR36.814
    PLdB = hPositioningPathLoss(radius(i), fc);
    PL = 10^(PLdB/10);
    
    % Add delay, pad and attenuate
    if ~isempty(rician_param)
        rx_tmp = [zeros(sampleDelay(i), 1); faded_signal(:, i); zeros(max(sampleDelay)-sampleDelay(i), 1)] ...
            /sqrt(PL);
    else
        rx_tmp = [zeros(sampleDelay(i), 1); tx_signal; zeros(max(sampleDelay)-sampleDelay(i), 1)] ...
            /sqrt(PL);
    end
    
    % awgn
    rx_signal(:, i) = awgn(rx_tmp, snr_db, 'measured', 'db');
    
end
size(rx_signal); % dimension = sample_length x sensor_length

% Plot received waveforms
if plot_signal
    plot_rx_signal(rx_signal, fs, snr_db);
end

%% Estimate Arrival Times
% ########## geolocation server do below procedures:
% (1) get received signal from each rf sensor, which may require large network bandwidth
%     it is assumed that rx cell array have already the received signal from all sensors 
% (2) choose reference sensor which is nearest to target emitter, 
%     so receieved signal of the sensor is largest, which will be used as reference signal
% (3) compute cross correlation between largest received signal and other received signal
% (4) estimate delay which is at point of maximum correlation

% choose reference sensor which is nearest to target
% ref_sensor: reference sensor number
[ref_sensor_number] = choose_reference_sensor(rx_signal);
ref_sensor_number
ref_signal = rx_signal(:, ref_sensor_number);

delayEst = zeros(1, sensor_length);
corr = zeros(nfft, sensor_length);
for i = 1 : sensor_length

    % Correlate received signal with reference sensor signal
    % ####### xcorr option = 'none'(default), 'biased', 'unbiased', 'coeff'
    % ####### xcorr option is not important, 
    % ####### but 'coeff' option is optimal?
    % ####### 'coeff' option: Normalizes the sequence so the autocorrelations at zero lag are identically 1.0.
%     c = abs(xcorr(rx_signal(:, i), ref_signal, 'unbiased')); 
%     c = abs(xcorr(rx_signal(:, i), ref_signal, 'coeff'));
    c = abs(xcorr(rx_signal(:, i), ref_signal));
    
    % Reduced length of correlation vector for positioning and plotting
    c(1:length(ref_signal)) = [];    % Remove meaningless result at beginning
    corr(:, i) = c(1 : nfft);   % Extract an OFDM symbol's worth of data
    
    % Delay estimate is at point of maximum correlation
    delayEst(i) = find(corr(:, i) == max(corr(:, i)));

end
size(corr);
corr;

% Plot correlation
if plot_signal
    plot_xcorrelation(corr, fs, ref_sensor_number);
end

%% Compute TDOA between reference sensor and other sensor, and plot constant TDOA hyperbolas

delayEst
delta_distance_meter = delayEst / fs * speedOfLight

if plot_position
    figure(H);
    for n = 1 : sensor_length
        % skip reference sensor itself
        if n == ref_sensor_number
            continue;
        end
        
        dd = (delayEst(n) - 1) / fs * speedOfLight;
        [x, y] = hPositioningTDOACurve(sensor_position(n, :), ...
            sensor_position(ref_sensor_number, :), dd);

        plot(x, y, 'k:', 'LineWidth', 2);
    end
    
    % Plot rf senosr and target emitter marker
    plot(target_position(1), target_position(2), ...
        'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'r', 'LineWidth', 2);
    
    grid on;

end

%%

% #### target_position must be 2 x 1
target_position = target_position';

[position_error_hyperbolic, position_error_linear, position_error_torrieri, ...
    x_est_torrieri, cep, gdop, lambda1, lambda2, theta_rad] = ...
    locate_target_tdoa(sensor_position, ref_sensor_number, delta_distance_meter, rx_signal, fs, ...
    snr_db, use_only_torrieri_method, target_position, bw_mhz);

end


