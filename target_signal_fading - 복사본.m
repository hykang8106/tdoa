function [] = target_signal_fading(plot_signal)

sensor_position = [3000, 2000];
target_position = [317, -51];

% plot_signal = 1;
snr_db =5;

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

% Plot transmit waveforms from target emitter
if plot_signal
    plot_tx_signal(tx_signal, fs);
%     my_hPositioningPlotTx(tx_signal, info.SamplingRate);
end
% fprintf('sample rate = %f MHz, distance per sample = %f meter\n', ...
%     info.SamplingRate / 1e6, physconst('LightSpeed') / info.SamplingRate);

% speedOfLight = 299792458.0; % Speed of light in m/s
speedOfLight = physconst('LightSpeed'); % Speed of light in m/s

% sampleDelay = zeros(1, sensor_length);
% radius = zeros(1, sensor_length);
% radius = cell(1, sensor_length);
% for i = 1 : sensor_length
%    [~, radius(i)] = cart2pol(sensor_position(i, 1) - target_position(1), ...
%        sensor_position(i, 2) - target_position(2));
%    delay = radius(i) / speedOfLight;                 % Delay in seconds
% %    [~, radius{i}] = cart2pol(rf_sensor{i}.Position(1), rf_sensor{i}.Position(2));
% %    delay = radius{i}/speedOfLight;                  % Delay in seconds
%    sampleDelay(i) = round(delay * info.SamplingRate); % Delay in samples 
% end
% sampleDelay;

[~, radius] = cart2pol(sensor_position(1) - target_position(1), sensor_position(2) - target_position(2));
radius
delay = radius / speedOfLight;                 % Delay in seconds
delay
%    [~, radius{i}] = cart2pol(rf_sensor{i}.Position(1), rf_sensor{i}.Position(2));
%    delay = radius{i}/speedOfLight;                  % Delay in seconds
sampleDelay = round(delay * info.SamplingRate); % Delay in samples
sampleDelay

% for i = 1 : sensor_length
%     % Urban Macro LOS path loss per TR36.814
%     PLdB = hPositioningPathLoss(radius(i), 2.1e9);
%     PL = 10^(PLdB/10);
%     
%     % Add delay, pad and attenuate
%     rx_tmp = [zeros(sampleDelay(i), 1); tx_signal; zeros(max(sampleDelay)-sampleDelay(i), 1)] ...
%         /sqrt(PL);
%     
%     % awgn
%     rx_signal(:, i) = awgn(rx_tmp, snr_db, 'measured', 'db');
% end

% Urban Macro LOS path loss per TR36.814
PLdB = hPositioningPathLoss(radius, 2.1e9)
PL = 10^(PLdB/10);

% Add delay, pad and attenuate
rx_tmp = [zeros(sampleDelay, 1); tx_signal; zeros(max(sampleDelay)-sampleDelay, 1)] ...
    /sqrt(PL);

% awgn
rx_signal = awgn(rx_tmp, snr_db, 'measured', 'db');

size(rx_signal); % dimension = sample_length x sensor_length

% % Plot received waveforms
% if plot_signal
%     plot_rx_signal(rx_signal, fs, snr_db);
% %     my_hPositioningPlotRx(rx, info.SamplingRate, snr_db);
% end

% ###################################################################################################
% ### Setting the maximum path delay greater than 100 samples may generate an ¡®Out of memory' error.
% ### when fs = 3.84e6 hz(sample period = 0.26 us), max path delay MUST BE LESS than 26 us
% ###################################################################################################

ts = 1 / fs;
fd = 0; % max doppler shift
% ##### rician k: ratio between specular power and diffuse power for a direct los path, 
% k = 1 ~ 10, rayleigh fading: k = 0, default: k = 1
% k = 1;
% k = 1;
k = [3, 1];

% ########## For outdoor environments, 
% ########## path delays after the first are typically between 100 ns(= 1e-7 s) and 10 us(= 1e-5 s)
% ########## 100 ns = 30 m, 10 us = 3 km
tau = [0, 7e-6];
% pdb = [-20, -30];
pdb = [0, -3];
chan = ricianchan(ts, fd, k, tau, pdb);
chan
% chan.PathDelays = 0;
% % chan.PathDelays = delay;
% chan.AvgPathGaindB = 0;
% chan.AvgPathGaindB = -20;
% chan.AvgPathGaindB = -PLdB;
% chan.NormalizePathGains = 0;
% ###############################################################################
% ##### StoreHistory must be false if MaxDopplerShift is zero.
% ##### channel visualization tool CANNOT BE USED when MaxDopplerShift is zero.
% ###############################################################################
% chan.StoreHistory = true; 
% chan

% % Ts = 1e-4;
% % fd = 100;
% k = [3 1];
% tau = [0 1e-5 1.5e-5 3e-5];
% pdb = [0 -1 -2 -2.5];
% h = ricianchan(Ts, fd, k, tau, pdb);

faded_rx_signal = filter(chan, tx_signal);
chan
size(faded_rx_signal)
% plot(chan);

if plot_signal
    % ###### snr_db input is dummy
    plot_rx_signal([tx_signal, faded_rx_signal], fs, snr_db);
%     plot_rx_signal(faded_rx_signal, fs, snr_db);
end

signal_compare = [tx_signal, faded_rx_signal];

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

end
