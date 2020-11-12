function [position_error_hyperbolic, position_error_linear, position_error_torrieri] = ...
    simulate_tdoa_target_fixed(rf_sensor_length, snr_db, randomize_sensor_distance, use_only_torrieri_method, ...
    plot_position, plot_signal)
% simulate target location using tdoa method
% 
% ################################################################
% ##### target emitter is fixed at (0,0)
% ##### sensors are randomly located within 4 ~ 10 km from target
% ################################################################
% 
% [input]
% - rf_sensor_length: number of rf sensor, 3 ~ 7
%   ofcom experience show sensor number greater than 5 dont improve location accuracy
% - snr_db: snr in db
% - randomize_sensor_distance: boolean. 0 = sensor 1 is nearest target, 1 = any sensor can be nearest target
% - use_only_torrieri_method:
%   when 0, use all method: torrieri(least square estimator), linear equation solver, hyperbolic equation solver
%   when 1, use only torrieri method
% - plot_position: boolean. to control plot of target, sensor, tdoa curve
% - plot_signal: boolean. to control plot of tx, rx, correlation
% [usage]
% simulate_tdoa_target_fixed(4, 5, 1, 0, 1, 1);

% #####################################################################
% ### CAUTION: 
% ### when writing code, dont forget target position be (0,0) 
% #####################################################################

%%
if rf_sensor_length > 7 || rf_sensor_length < 3
    fprintf(2, '##### sensor number: 3 ~ 7\n');
    return;
end

rf_sensor = cell(1,rf_sensor_length);

%% Transmitter(target emitter: lte base station) Configuration

% ### must comment out in function implementation, only useful in script implementation
% ### run learn_rng.m
% rng('default'); % Initialize the random number generator stream

% target emitter to be located is lte base station
enb = lteRMCDL('R.5');   % Get configuration based on RMC
enb.NCellID = 0;     % Set unique cell identity
enb.TotSubframes = 1;  % Number of subframes to generate
enb.NPRSRB = 2;        % Bandwidth of PRS in resource blocks
enb.IPRS = 0;          % PRS configuration index
enb.PRSPeriod = 'On';  % PRS present in all subframes
enb;

% get rf sensor position, randomly located
for i=1:rf_sensor_length
    rf_sensor{i}.Position = hPositioningPosition(i-1, rf_sensor_length);
    rf_sensor{i}.Position;
end

% randomize sensor distance from target
% in original version(hPositioningPosition.m), 
% sensor 1 was nearest from target, sensor 4 was furthest from target
if randomize_sensor_distance
    rf_sensor = randomize_sensor_distance_from_target(rf_sensor);
end

% get ofdm modulation related information
info = lteOFDMInfo(enb);

%% Plot Location of rf sensor and target emitter  
if plot_position
    H = my_hPositioningPlotPositions(rf_sensor);
end

%% Generate Transmissions from target emitter and Plot Transmitted Waveforms

rb_grid = lteDLResourceGrid(enb);          % Empty resource grid
rb_grid(ltePRSIndices(enb)) = ltePRS(enb); % Map PRS REs
tx{1} = lteOFDMModulate(enb, rb_grid);        % OFDM modulate

% Plot transmit waveforms from target emitter
if plot_signal
    my_hPositioningPlotTx(tx, info.SamplingRate);
end
% fprintf('sample rate = %f MHz\n', info.SamplingRate / 1e6);

%% Compute delays from target emitter to rf sensors
speedOfLight = 299792458.0; % Speed of light in m/s

sampleDelay = zeros(1, rf_sensor_length);
radius = cell(1, rf_sensor_length);
for i = 1:rf_sensor_length
   [~, radius{i}] = cart2pol(rf_sensor{i}.Position(1), rf_sensor{i}.Position(2));
   delay = radius{i}/speedOfLight;                  % Delay in seconds
   sampleDelay(i) = round(delay*info.SamplingRate); % Delay in samples 
end

%% Create Received Waveforms in each rf sensor and Plot Received Waveform

rx = cell(1, rf_sensor_length);
for i = 1:rf_sensor_length
    % Urban Macro LOS path loss per TR36.814
    PLdB = hPositioningPathLoss(radius{i}, 2.1e9);
    PL = 10^(PLdB/10);
    
    % Add delay, pad and attenuate
    rx{i} = [zeros(sampleDelay(i), 1); tx{1}; zeros(max(sampleDelay)-sampleDelay(i), 1)] ...
        /sqrt(PL);
    
    % awgn
    rx{i} = awgn(rx{i}, snr_db, 'measured', 'db');
end

% Plot received waveforms
if plot_signal
    my_hPositioningPlotRx(rx, info.SamplingRate, snr_db);
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
[ref_sensor_number] = choose_reference_sensor(rx);
ref_signal = rx{ref_sensor_number};

corr = cell(1, rf_sensor_length);
delayEst = zeros(1, rf_sensor_length);
% corr = cell(1, rf_sensor_length - 1);
% delayEst = zeros(1, rf_sensor_length - 1);

for i = 1:rf_sensor_length
    % #### auto-correlation of reference rf sensor is NOT necessary
    % #### but this is included for simple code, and later use(auto-correlation is equal to fft)
%     if i == ref_sensor
%         continue;   % skip referenec sensor
%     end

    % Correlate received signal with reference sensor signal
    c = abs(xcorr(rx{i},ref_signal));
    
    % Reduced length of correlation vector for positioning and plotting
    c(1:length(ref_signal)) = [];    % Remove meaningless result at beginning
    corr{i} = c(1:info.Nfft);   % Extract an OFDM symbol's worth of data
    
    % Delay estimate is at point of maximum correlation
    delayEst(i) = find(corr{i} == max(corr{i}));
end

% Plot correlation
if plot_signal
    my_hPositioningPlotCorr(corr, info.SamplingRate, ref_sensor_number);
end

%% Compute TDOA between reference sensor and other sensor, and plot constant TDOA hyperbolas

delayEst
delta_distance_meter = delayEst / info.SamplingRate * speedOfLight

if plot_position
    figure(H);
    for n = 1 : rf_sensor_length
        % skip reference sensor itself
        if n == ref_sensor_number
            continue;
        end
        
        dd = (delayEst(n) - 1) / info.SamplingRate * speedOfLight;
        [x, y] = hPositioningTDOACurve(rf_sensor{n}.Position, ...
            rf_sensor{ref_sensor_number}.Position, dd);
        plot(x, y, 'k:', 'LineWidth', 2);
    end
    
    % Plot rf senosr and target emitter marker
    plot(0, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'r', 'LineWidth', 2);
    % plot(0, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'c', 'LineWidth', 2);
    
    grid on;
end

%%

% #### rf_sensor is cell array
% #### convert to matrix for sensor_position input of locate_target_tdoa function

% convert cell to sensor position matrix(dimension = sensor_length x 2)
sensor_length = length(rf_sensor);
sensor_position = zeros(sensor_length, 2);
for n = 1 : sensor_length
    rf_sensor{n}.Position;
    sensor_position(n, :) = rf_sensor{n}.Position;
end
sensor_position

% ##### must be 2 x 1
target_position = [0; 0];

% ##############################################################################
% #### ndlrb = 15, nprsrb = 2, subframe_length = 1
% #### [bw_mhz, fs, nfft, sample_length] = get_bw_from_prs_spec_db(15, 2, 1)
% ##############################################################################
bw_mhz = 0.4;
[position_error_hyperbolic, position_error_linear, position_error_torrieri, ...
    x_est_torrieri, cep, gdop, lambda1, lambda2, theta_rad] = ...
    locate_target_tdoa(sensor_position, ref_sensor_number, delta_distance_meter, rx, info.SamplingRate, snr_db, ...
    use_only_torrieri_method, target_position, bw_mhz);

if plot_position
    % make title text
    title_base = sprintf('[location error, SNR %d dB]', snr_db);
    if use_only_torrieri_method
        title_text = sprintf('%s (LS) %d m', title_base, round(position_error_torrieri));
    else
        if sensor_length == 3
            title_text = sprintf('%s (LS) %d m, (HE) %d m', ...
                title_base, round(position_error_torrieri), round(position_error_hyperbolic));
        else
            title_text = sprintf('%s (LS) %d m, (LE) %d m', ...
                title_base, round(position_error_torrieri), round(position_error_linear));
        end
    end
    % show title
    title(gca, title_text, 'color', 'b', 'fontweight', 'bold');
end

end

