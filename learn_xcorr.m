function [] = learn_xcorr(plot_signal)

snr_db = 5;

sensor_position = ...
    [1.5000000e+03   3.5000000e+03
    -1.5000000e+03   2.0000000e+03
    -2.0000000e+03  -2.0000000e+03
     2.0000000e+03  -2.5000000e+03];
 
target_position = [0, 0];

sensor_length = size(sensor_position, 1);

% target emitter to be located is lte base station
enb = lteRMCDL('R.5');   % Get configuration based on RMC
enb.NCellID = 0;     % Set unique cell identity
enb.TotSubframes = 1;  % Number of subframes to generate
enb.NPRSRB = 2;        % Bandwidth of PRS in resource blocks
enb.IPRS = 0;          % PRS configuration index
enb.PRSPeriod = 'On';  % PRS present in all subframes
enb

% get ofdm modulation related information
info = lteOFDMInfo(enb);

resource_grid = lteDLResourceGrid(enb);          % Empty resource grid
resource_grid(ltePRSIndices(enb)) = ltePRS(enb); % Map PRS REs
tx_signal = lteOFDMModulate(enb, resource_grid);        % OFDM modulate
size(tx_signal) % column vector

% Plot transmit waveforms from target emitter
if plot_signal
    plot_tx_signal(tx_signal, info.SamplingRate);
%     my_hPositioningPlotTx(tx_signal, info.SamplingRate);
end

% fprintf('sample rate = %f MHz, distance per sample = %f meter\n', ...
%     info.SamplingRate / 1e6, physconst('LightSpeed') / info.SamplingRate);

speedOfLight = physconst('LightSpeed'); % Speed of light in m/s

sampleDelay = zeros(1, sensor_length);
radius = zeros(1, sensor_length);
% radius = cell(1, sensor_length);
for i = 1 : sensor_length
   [~, radius(i)] = cart2pol(sensor_position(i, 1) - target_position(1), ...
       sensor_position(i, 2) - target_position(2));
   delay = radius(i) / speedOfLight;                 % Delay in seconds
%    [~, radius{i}] = cart2pol(rf_sensor{i}.Position(1), rf_sensor{i}.Position(2));
%    delay = radius{i}/speedOfLight;                  % Delay in seconds
   sampleDelay(i) = round(delay * info.SamplingRate); % Delay in samples 
end
sampleDelay

for i = 1 : sensor_length
    % Urban Macro LOS path loss per TR36.814
    PLdB = hPositioningPathLoss(radius(i), 2.1e9);
    PL = 10^(PLdB/10);
    
    % Add delay, pad and attenuate
    rx_tmp = [zeros(sampleDelay(i), 1); tx_signal; zeros(max(sampleDelay)-sampleDelay(i), 1)] ...
        /sqrt(PL);
    
    % awgn
    rx_signal(:, i) = awgn(rx_tmp, snr_db, 'measured', 'db');
end
size(rx_signal) % dimension = sample_length x sensor_length

if plot_signal
    plot_rx_signal(rx_signal, info.SamplingRate, snr_db);
end

sensor_pair = nchoosek(1 : sensor_length, 2);
pair_length = size(sensor_pair, 1);
delayEst = zeros(1, pair_length);

corr = zeros(2 * size(rx_signal, 1) - 1, pair_length);

for i = 1 : pair_length 
    sp = sensor_pair(i, :);
    plus_sensor_number = sp(1);
    minus_sensor_number = sp(2);
    
    % Correlate received signal with reference sensor signal
    corr(:, i) = abs(xcorr(rx_signal(:, plus_sensor_number), rx_signal(:, minus_sensor_number)));
    
%     % Reduced length of correlation vector for positioning and plotting
%     c(1 : size(rx_signal, 1)) = [];    % Remove meaningless result at beginning
%     corr(:, i) = c(1 : info.Nfft);   % Extract an OFDM symbol's worth of data
%     
%     % Delay estimate is at point of maximum correlation
%     delayEst(i) = find(corr(:, i) == max(corr(:, i)));
end
size(corr)

[Y, I] = max(corr)

delayEst = I - size(rx_signal, 1)

figure;
plot(corr, '.-');
grid on;

% for i = 1 : pair_length 
%     sp = sensor_pair(i, :);
%     plus_sensor_number = sp(1);
%     minus_sensor_number = sp(2);
%     
%     % Correlate received signal with reference sensor signal
%     c = abs(xcorr(rx_signal(:, plus_sensor_number), rx_signal(:, minus_sensor_number)));
%     
%     % Reduced length of correlation vector for positioning and plotting
%     c(1 : size(rx_signal, 1)) = [];    % Remove meaningless result at beginning
%     corr(:, i) = c(1 : info.Nfft);   % Extract an OFDM symbol's worth of data
%     
%     % Delay estimate is at point of maximum correlation
%     delayEst(i) = find(corr(:, i) == max(corr(:, i)));
% end
% size(corr);
% 
% delayEst
% delta_distance_meter = delayEst / info.SamplingRate * speedOfLight

end