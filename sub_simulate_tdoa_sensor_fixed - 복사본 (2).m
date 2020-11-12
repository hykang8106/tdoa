function [position_error_torrieri, x_est_torrieri, cep, gdop, lambda1, lambda2, theta_rad] = ...
    sub_simulate_tdoa_sensor_fixed(sensor_position, snr_db, target_position, ...
    plot_signal, plot_position, use_only_torrieri_method, target_radius_meter, rician_param, ...
    target_signal_spec)
% simulate target location using tdoa method
%
% this differ from simulate_tdoa_target_fixed:
% ######################################################################
% ##### target position is given by input or randomly determined
% ##### sensors location are fixed: uca + (0,0)
% ##### only torrieri algorithm(least square estimator) is implemented
% ##### array variable is used instead of cell variable
% ######################################################################
%
% [input]
% - sensor_length: number of sensor, 3 ~ 7
% - sensor_exist_at_uca_center: boolean
% - uca_radius_meter: uca radius in meter
% - snr_db: snr in db
% - target_position: target position. if empty, random position
% - radius_ratio: 
%   ### when radius_ratio <= 1, target is located inside uca
%   ### when radius_ratio > 1, target can be located outside uca
% - sensor_position_filename: if not empty, sensor position in file is used, ascii file, .txt
%   ### MUST USE [], not implemented yet
% [usage]
% simulate_tdoa_sensor_fixed(4, 1, 7e3, 5, [], 1, []);
% simulate_tdoa_sensor_fixed(4, 1, 7e3, 5, [1055, -3000], 1, []);
%

% #########################################################################
% ### CAUTION: 
% ### when writing code, dont forget target position NOT be (0,0) 
% #########################################################################

%% Check input parameter

sensor_length = size(sensor_position, 1);

% if sensor_length > 7 || sensor_length < 3
%     fprintf('sensor number: 3 ~ 7\n');
%     return;
% end
% % rf_sensor = cell(1,sensor_length);
% 
% % sensor_exist_at_uca_center = 1;
% 
% % when 0, use all method: torrieri(least square estimator), linear equation solver, hyperbolic equation solver
% % when 1, use only torrieri method
% use_only_torrieri_method = 1;
% 
% % radius_ratio = 1.0; % when radius_ratio <= 1.0, target is located inside uca
% % radius_ratio = 1.5;
% % radius_ratio = 2.0; % when radius_ratio > 1.0, target can be located outside uca
% 
% if ~isempty(target_position)
%     if target_position(1) > (uca_radius_meter * radius_ratio) || ...
%         target_position(1) < -(uca_radius_meter * radius_ratio) || ...
%         target_position(2) > (uca_radius_meter * radius_ratio) || ...
%         target_position(2) < -(uca_radius_meter * radius_ratio)
%         fprintf('###### target position is out of range\n');
%         return;
%     end
% end
% 
% plot_signal = 0; % control plot of tx, rx, correlation
% plot_position = 0; % control plot of target, sensor, tdoa curve

%% Plot Location of rf sensor and target emitter  
if plot_position
    H = plot_sensor_position(sensor_position, target_position, target_radius_meter);
end

%% Transmitter(target emitter: lte base station) Configuration

ndlrb = target_signal_spec(1);
nprsrb = target_signal_spec(2);
subframe_length = target_signal_spec(3);

[bw_mhz, fs, nfft, sample_length] = get_bw_from_prs_spec_db(ndlrb, nprsrb, subframe_length);

[tx_signal, fs, nfft] = generate_target_signal_lte_prs(ndlrb, nprsrb, subframe_length);
tx_signal_length = length(tx_signal);

% ############################## 161029 ############################################################

% % ### must comment out in function implementation, only useful in script implementation
% % ### run learn_rng.m
% % rng('default'); % Initialize the random number generator stream
% 
% % target emitter to be located is lte base station
% enb = lteRMCDL('R.5');   % Get configuration based on RMC
% enb.NCellID = 0;     % Set unique cell identity
% enb.TotSubframes = 1;  % Number of subframes to generate
% enb.NPRSRB = 2;        % Bandwidth of PRS in resource blocks
% enb.IPRS = 0;          % PRS configuration index
% enb.PRSPeriod = 'On';  % PRS present in all subframes
% 
% % % sensor_position, dimension = sensor_length x 2
% % % #################################################################################
% % % ### write code to get sensor position from file which have sensor position
% % % ### this is good for sensor planning in real environment
% % % #################################################################################
% % if ~isempty(sensor_position_filename)
% %     sensor_position = load(sensor_position_filename);
% % else
% %     sensor_position = get_sensor_position(sensor_length, uca_radius_meter, sensor_exist_at_uca_center);
% % end
% % 
% % % % get rf sensor position, randomly located
% % % for i=1:sensor_length
% % %     rf_sensor{i}.Position = hPositioningPosition(i-1, sensor_length);
% % %     rf_sensor{i}.Position;
% % % end
% % 
% % if isempty(target_position)
% %     % target_position, dimension = 1 x 2
% %     target_position = random_target_position(uca_radius_meter, radius_ratio);
% % end
% 
% % % randomize sensor distance from target
% % % in original version(hPositioningPosition.m), 
% % % sensor 1 was nearest from target, sensor 4 was furthest from target
% % if randomize_sensor_distance
% %     rf_sensor = randomize_sensor_distance_from_target(rf_sensor);
% % end
% 
% % get ofdm modulation related information
% info = lteOFDMInfo(enb);
% 
% %% Plot Location of rf sensor and target emitter  
% % if plot_position
% %     H = plot_sensor_position(sensor_position, target_position, target_radius_meter);
% % end
% 
% %% Generate Transmissions from target emitter and Plot Transmitted Waveforms
% 
% % tx = cell(1,1);
% resource_grid = lteDLResourceGrid(enb);          % Empty resource grid
% resource_grid(ltePRSIndices(enb)) = ltePRS(enb); % Map PRS REs
% tx_signal = lteOFDMModulate(enb, resource_grid);        % OFDM modulate
% size(tx_signal); % column vector
% tx_signal_length = length(tx_signal);
% ################################# 161029 ##############################################################

% Plot transmit waveforms from target emitter
if plot_signal
    plot_tx_signal(tx_signal, fs);
%     plot_tx_signal(tx_signal, info.SamplingRate);
%     my_hPositioningPlotTx(tx_signal, info.SamplingRate);
end
fprintf('sample rate = %f MHz, distance per sample = %f meter\n', ...
    fs / 1e6, physconst('LightSpeed') / fs);

%% Compute delays from target emitter to rf sensors
% speedOfLight = 299792458.0; % Speed of light in m/s
speedOfLight = physconst('LightSpeed'); % Speed of light in m/s

sampleDelay = zeros(1, sensor_length);
radius = zeros(1, sensor_length);
% radius = cell(1, sensor_length);
faded_signal = zeros(tx_signal_length, sensor_length);
for i = 1 : sensor_length
    [~, radius(i)] = cart2pol(sensor_position(i, 1) - target_position(1), ...
        sensor_position(i, 2) - target_position(2));
    delay = radius(i) / speedOfLight;                 % Delay in seconds
    %    [~, radius{i}] = cart2pol(rf_sensor{i}.Position(1), rf_sensor{i}.Position(2));
    %    delay = radius{i}/speedOfLight;                  % Delay in seconds
    
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
    
    sampleDelay(i) = round(delay * fs); % Delay in samples
end
sampleDelay;

%% Create Received Waveforms in each rf sensor and Plot Received Waveform

% rx = cell(1, sensor_length);
rx_signal = zeros(tx_signal_length + max(sampleDelay), sensor_length);
for i = 1 : sensor_length
    % Urban Macro LOS path loss per TR36.814
    PLdB = hPositioningPathLoss(radius(i), 2.1e9);
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
    
    %     i = i + 1;
    
    %     % Add delay, pad and attenuate
    %     rx{i} = [zeros(sampleDelay(i), 1); tx_signal{1}; zeros(max(sampleDelay)-sampleDelay(i), 1)] ...
    %         /sqrt(PL);
    %
    %     % awgn
    %     rx{i} = awgn(rx{i}, snr_db, 'measured', 'db');
end
size(rx_signal); % dimension = sample_length x sensor_length

% Plot received waveforms
if plot_signal
    plot_rx_signal(rx_signal, fs, snr_db);
%     my_hPositioningPlotRx(rx, info.SamplingRate, snr_db);
end

% figure;
% plot(10*log10(abs(fft(rx{1}))));

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
ref_signal = rx_signal(:, ref_sensor_number);

% corr = cell(1, sensor_length);
delayEst = zeros(1, sensor_length);
% corr = cell(1, rf_sensor_length - 1);
% delayEst = zeros(1, rf_sensor_length - 1);
for i = 1 : sensor_length
    % #### auto-correlation of reference rf sensor is NOT necessary
    % #### but this is included for simple code, and later use(auto-correlation is equal to fft)
%     if i == ref_sensor
%         continue;   % skip referenec sensor
%     end

    % Correlate received signal with reference sensor signal
    c = abs(xcorr(rx_signal(:, i), ref_signal));
%     c = abs(xcorr(rx_signal{i},ref_signal));
    
    % Reduced length of correlation vector for positioning and plotting
    c(1:length(ref_signal)) = [];    % Remove meaningless result at beginning
    corr(:, i) = c(1 : nfft);   % Extract an OFDM symbol's worth of data
%     corr{i} = c(1:info.Nfft);   % Extract an OFDM symbol's worth of data
    
    % Delay estimate is at point of maximum correlation
    delayEst(i) = find(corr(:, i) == max(corr(:, i)));
%     delayEst(i) = find(corr{i} == max(corr{i}));
end
size(corr);

% Plot correlation
if plot_signal
    plot_xcorrelation(corr, fs, ref_sensor_number);
%     my_hPositioningPlotCorr(corr, info.SamplingRate, ref_sensor_number);
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
        
%         [x, y] = hPositioningTDOACurve(rf_sensor{n}.Position, ...
%             rf_sensor{ref_sensor_number}.Position, dd);

        plot(x, y, 'k:', 'LineWidth', 2);
    end
    
    % % Estimate time difference of arrival from each sensor
    % tdoa = hPositioningTDOA(delayEst,info.SamplingRate);
    %
    % % Plot hyperbolas
    % figure(1);
    % for j = 1:1
    %     for i = (j+1):rf_sensor_length
    %         dd = tdoa(i,j)*speedOfLight; % Delay distance
    %         [x, y] = hPositioningTDOACurve(rf_sensor{i}.Position, ...
    %             rf_sensor{j}.Position, dd);
    %         plot(x, y, 'k:', 'LineWidth', 2);
    %     end
    % end
    
    % for j = 1:rf_sensor_length
    %     for i = (j+1):rf_sensor_length
    %         dd = tdoa(i,j)*speedOfLight; % Delay distance
    %         [x, y] = hPositioningTDOACurve(rf_sensor{i}.Position, ...
    %             rf_sensor{j}.Position, dd);
    %         plot(x, y, 'k:', 'LineWidth', 2);
    %     end
    % end
    
    % Plot rf senosr and target emitter marker
    plot(target_position(1), target_position(2), ...
        'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'r', 'LineWidth', 2);
    
    grid on;
%     plot(0, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'r', 'LineWidth', 2);
    % plot(0, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'c', 'LineWidth', 2);
end

%%

% #### must be 2 x 1
target_position = target_position';

[position_error_hyperbolic, position_error_linear, position_error_torrieri, ...
    x_est_torrieri, cep, gdop, lambda1, lambda2, theta_rad] = ...
    locate_target_tdoa(sensor_position, ref_sensor_number, delta_distance_meter, rx_signal, fs, ...
    snr_db, use_only_torrieri_method, target_position, bw_mhz);

end


