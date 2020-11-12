function [] = target_signal_fading(plot_signal, signal_fading, snr_db)
% learn target signal fading
%
% [usage]
% target_signal_fading(1, 1, 0)
%

sensor_position_filename = 'sensor_position.txt';
rician_param_filename = 'rician_param.xlsx';

[sensor_position, uca_radius_meter] = get_sensor_position_from_file(sensor_position_filename);
sensor_length = size(sensor_position, 1);

% sensor_position = [3000, 2000];
target_position = [317, -51];

xlsrange = 'b3:d7'; % see excel file(blue bold cell)
rician_param = get_rician_parameter_from_excel_file(rician_param_filename, xlsrange);
if size(rician_param, 1) ~= sensor_length
    fprintf('##### row length of rician parameter must be same as sensor length\n');
    return;
end

% plot_signal = 1;
% snr_db =5;

% target emitter to be located is lte base station
enb = lteRMCDL('R.5');   % Get configuration based on RMC
enb.NCellID = 0;     % Set unique cell identity
enb.TotSubframes = 1;  % Number of subframes to generate
enb.NPRSRB = 2;        % Bandwidth of PRS in resource blocks
enb.IPRS = 0;          % PRS configuration index
enb.PRSPeriod = 'On';  % PRS present in all subframes

% get ofdm modulation related information
info = lteOFDMInfo(enb);
fs = info.SamplingRate;

% tx = cell(1,1);
resource_grid = lteDLResourceGrid(enb);          % Empty resource grid
resource_grid(ltePRSIndices(enb)) = ltePRS(enb); % Map PRS REs
tx_signal = lteOFDMModulate(enb, resource_grid);        % OFDM modulate
size(tx_signal); % column vector
tx_signal_length = length(tx_signal);

% Plot transmit waveforms from target emitter
% if plot_signal
%     plot_tx_signal(tx_signal, fs);
% %     my_hPositioningPlotTx(tx_signal, info.SamplingRate);
% end
% fprintf('sample rate = %f MHz, distance per sample = %f meter\n', ...
%     info.SamplingRate / 1e6, physconst('LightSpeed') / info.SamplingRate);

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
    
    k = rician_param{i, 1};
    if signal_fading
        if ~isempty(k)
            % tau in rician_param is delta tau ratio(not absolute delay), see excel file
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
    
    sampleDelay(i) = round(delay * info.SamplingRate); % Delay in samples
end
sampleDelay;

% [~, radius] = cart2pol(sensor_position(1) - target_position(1), sensor_position(2) - target_position(2));
% radius
% delay = radius / speedOfLight;                 % Delay in seconds
% delay
% %    [~, radius{i}] = cart2pol(rf_sensor{i}.Position(1), rf_sensor{i}.Position(2));
% %    delay = radius{i}/speedOfLight;                  % Delay in seconds
% sampleDelay = round(delay * info.SamplingRate); % Delay in samples
% sampleDelay

rx_signal = zeros(tx_signal_length + max(sampleDelay), sensor_length);
for i = 1 : sensor_length
    % Urban Macro LOS path loss per TR36.814
    PLdB = hPositioningPathLoss(radius(i), 2.1e9);
    PL = 10^(PLdB/10);
    
    if signal_fading
        % Add delay, pad and attenuate
        rx_tmp = [zeros(sampleDelay(i), 1); faded_signal(:, i); zeros(max(sampleDelay)-sampleDelay(i), 1)] ...
            /sqrt(PL);
    else
        % Add delay, pad and attenuate
        rx_tmp = [zeros(sampleDelay(i), 1); tx_signal; zeros(max(sampleDelay)-sampleDelay(i), 1)] ...
            /sqrt(PL);
    end
    
    % awgn
    rx_signal(:, i) = awgn(rx_tmp, snr_db, 'measured', 'db');
end

% ###################################################################################################
% ### Setting the maximum path delay greater than 100 samples may generate an ¡®Out of memory' error.
% ### when fs = 3.84e6 hz(sample period = 0.26 us), max path delay MUST BE LESS than 26 us
% ###################################################################################################

% ########## For outdoor environments, 
% ########## path delays after the first are typically between 100 ns(= 1e-7 s) and 10 us(= 1e-5 s)
% ########## 100 ns = 30 m, 10 us = 3 km
% tau = [0, 7e-6];
% pdb = [-20, -30];
% pdb = [0, -3];

% ##### rician k: ratio between specular power and diffuse power for a direct los path, 
% k = 1 ~ 10, rayleigh fading: k = 0, default: k = 1
% k = 1;
% k = 1;
% k = [3, 1];

% % tau: path delay in second
% 1 direct los path, 2 reflected path(rayleigh)
% tau = [0, 1e-6, 2e-6];
% % pdb: average path gain in db
% pdb = [0, -10, -20];
% chan = ricianchan(ts,fd,k,tau,pdb);

% Ts = 1e-4;
% fd = 100;
% kFactor = [3, 1, 0.2, 0];
% tau = [0, 1e-5, 1.5e-5, 3e-5];
% pdb = [0, -1, -2, -2.5];
% h = ricianchan(Ts, fd, kFactor, tau, pdb);

% k = [3, 1];
% tau = [0, 1e-6];
% pdb = [0, -3];
% [faded_signal, chan] = rician_fading(tx_signal, fs, k, tau, pdb);
% 
% if plot_signal
%     % ###### snr_db input is dummy
%     plot_rx_signal([tx_signal, faded_signal], fs, snr_db);
% %     plot_rx_signal(faded_rx_signal, fs, snr_db);
% end
% 
% % Urban Macro LOS path loss per TR36.814
% PLdB = hPositioningPathLoss(radius, 2.1e9)
% PL = 10^(PLdB/10);
% 
% % Add delay, pad and attenuate
% rx_tmp = [zeros(sampleDelay, 1); faded_signal; zeros(max(sampleDelay)-sampleDelay, 1)] ...
%     /sqrt(PL);
% % rx_tmp = [zeros(sampleDelay, 1); tx_signal; zeros(max(sampleDelay)-sampleDelay, 1)] ...
% %     /sqrt(PL);
% 
% % awgn
% rx_signal = awgn(rx_tmp, snr_db, 'measured', 'db');
% 
% size(rx_signal); % dimension = sample_length x sensor_length

% Plot received waveforms
if plot_signal
    plot_rx_signal(rx_signal, fs, snr_db);
%     my_hPositioningPlotRx(rx, info.SamplingRate, snr_db);
end

% plot(chan);

% if plot_signal
%     % ###### snr_db input is dummy
%     plot_rx_signal([tx_signal, faded_rx_signal], fs, snr_db);
% %     plot_rx_signal(faded_rx_signal, fs, snr_db);
% end

% signal_compare = [tx_signal, faded_rx_signal];

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
    corr(:, i) = c(1 : info.Nfft);   % Extract an OFDM symbol's worth of data
%     corr{i} = c(1:info.Nfft);   % Extract an OFDM symbol's worth of data
    
    % Delay estimate is at point of maximum correlation
    delayEst(i) = find(corr(:, i) == max(corr(:, i)));
%     delayEst(i) = find(corr{i} == max(corr{i}));
end
size(corr);
delayEst

% Plot correlation
if plot_signal
    plot_xcorrelation(corr, info.SamplingRate, ref_sensor_number);
%     my_hPositioningPlotCorr(corr, info.SamplingRate, ref_sensor_number);
end

end

%%
% function [faded_signal, chan] = rician_fading(tx_signal, fs, k, tau, pdb)
% 
% % ###################################################################################################
% % ### Setting the maximum path delay greater than 100 samples may generate an ¡®Out of memory' error.
% % ### when fs = 3.84e6 hz(sample period = 0.26 us), max path delay MUST BE LESS than 26 us
% % ###################################################################################################
% 
% ts = 1 / fs;
% fd = 0; % max doppler shift
% % ##### rician k: ratio between specular power and diffuse power for a direct los path, 
% % k = 1 ~ 10, rayleigh fading: k = 0, default: k = 1
% % k = 1;
% % k = 1;
% % k = [3, 1];
% 
% % ########## For outdoor environments, 
% % ########## path delays after the first are typically between 100 ns(= 1e-7 s) and 10 us(= 1e-5 s)
% % ########## 100 ns = 30 m, 10 us = 3 km
% % tau = [0, 7e-6];
% % pdb = [-20, -30];
% % pdb = [0, -3];
% chan = ricianchan(ts, fd, k, tau, pdb);
% chan;
% % chan.PathDelays = 0;
% % % chan.PathDelays = delay;
% % chan.AvgPathGaindB = 0;
% % chan.AvgPathGaindB = -20;
% % chan.AvgPathGaindB = -PLdB;
% % chan.NormalizePathGains = 0;
% % ###############################################################################
% % ##### StoreHistory must be false if MaxDopplerShift is zero.
% % ##### channel visualization tool CANNOT BE USED when MaxDopplerShift is zero.
% % ###############################################################################
% % chan.StoreHistory = true; 
% % chan
% 
% % % Ts = 1e-4;
% % % fd = 100;
% % k = [3 1];
% % tau = [0 1e-5 1.5e-5 3e-5];
% % pdb = [0 -1 -2 -2.5];
% % h = ricianchan(Ts, fd, k, tau, pdb);
% 
% faded_signal = filter(chan, tx_signal);
% chan;
% size(faded_signal);
% 
% if chan.ChannelFilterDelay
%     faded_signal = [faded_signal(chan.ChannelFilterDelay + 1 : end); zeros(chan.ChannelFilterDelay, 1)];
% end
% 
% end


