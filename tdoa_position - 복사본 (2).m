function [] = tdoa_position(rf_sensor_length, snr_db)
% target location using tdoa method
% modified to simulate tdoa location of target emitter used for spectrum monitoring
% original was PRSPositioningExample.m in LTE System Toolbox example
% target emitter is located at (0,0)
% 
% [input]
% - rf_sensor_length: number of rf sensor, 3 ~ 5
%   ofcom experience show sensor number greater than 5 dont improve location accuracy
% - snr_db: snr in db
% [usage]
% tdoa_position(4, 0)
%

%%
% ######################## begin of my code #############

% default
randomize_sensor_distance = 1;

if rf_sensor_length > 5 || rf_sensor_length < 3
    fprintf('rf sensor number: 3 ~ 5\n');
    return;
end
rf_sensor = cell(1,rf_sensor_length);

%% Transmitter(target emitter: lte base station) Configuration

% ### must comment out in function implementation, useful in script implementation
% ### run learn_rng.m
% rng('default'); % Initialize the random number generator stream

% target emitter to be located is lte base station
enb = lteRMCDL('R.5');   % Get configuration based on RMC
enb.NCellID = 0;     % Set unique cell identity
enb.TotSubframes = 1;  % Number of subframes to generate
enb.NPRSRB = 2;        % Bandwidth of PRS in resource blocks
enb.IPRS = 0;          % PRS configuration index
enb.PRSPeriod = 'On';  % PRS present in all subframes

% get rf sensor position, randomly located
for i=1:rf_sensor_length
    rf_sensor{i}.Position = hPositioningPosition(i-1, rf_sensor_length);
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

H = my_hPositioningPlotPositions(rf_sensor);

%% Generate Transmissions from target emitter and Plot Transmitted Waveforms

% tx = cell(1,1);
grid = lteDLResourceGrid(enb);          % Empty resource grid
grid(ltePRSIndices(enb)) = ltePRS(enb); % Map PRS REs
tx{1} = lteOFDMModulate(enb, grid);        % OFDM modulate

% Plot transmit waveforms from target emitter
my_hPositioningPlotTx(tx, info.SamplingRate);

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
my_hPositioningPlotRx(rx, info.SamplingRate);

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
[ref_sensor] = choose_reference_sensor(rx);
ref_signal = rx{ref_sensor};

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
my_hPositioningPlotCorr(corr, info.SamplingRate, ref_sensor);

%% Compute TDOA between reference sensor and other sensor, and plot constant TDOA hyperbolas

delayEst

figure(H);
for n = 1 : rf_sensor_length
    % skip reference sensor itself
    if n == ref_sensor
        continue;
    end
    
    dd = (delayEst(n) - 1) / info.SamplingRate * speedOfLight;
    [x, y] = hPositioningTDOACurve(rf_sensor{n}.Position, ...
        rf_sensor{ref_sensor}.Position, dd);
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
plot(0, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'c', 'LineWidth', 2);

%%
return;

% ######################### end of my code ###################

%% Time Difference Of Arrival Positioning Using PRS
% This example shows how to use the Time Difference Of Arrival (TDOA)
% positioning approach in conjunction with the Release 9 Positioning
% Reference Signal (PRS) to calculate the position of a User Equipment (UE)
% within a network of eNodeBs using the LTE System Toolbox(TM).

% Copyright 2011-2013 The MathWorks, Inc.

%% Introduction
% In this example, a number of eNodeB transmissions are created and
% combined with different delays and received powers to model the reception
% of all the eNodeB waveforms by one UE. The UE performs correlation with
% the Positioning Reference Signal (PRS) to establish the delay from each
% eNodeB and subsequently the delay difference between all pairs of
% eNodeBs. These delay differences are used to compute hyperbolas of
% constant delay difference, which are plotted relative to known eNodeB
% positions and intersect at the position of the UE.

%% Transmitter Configuration
% A set (cell array) of eNodeB configurations |enb| is created, with the
% number of eNodeBs specified by |NeNodeB|. The configurations are derived
% from Reference Measurement Channel (RMC) R.5 using
% <matlab:doc('lteRMCDL') lteRMCDL>. R.5 describes a 3 MHz bandwidth
% DownLink Shared Channel (PDSCH) transmission using 64-QAM modulation. For
% each eNodeB the configuration is updated to make the cell identity
% |NCellID| unique and the PRS parameters |NPRSRB|, |IPRS| and |PRSPeriod|
% are set.
%
% A random position given by an X and Y coordinate for each eNodeB is
% generated by <matlab:edit('hPositioningPosition.m')
% hPositioningPosition.m>.

rng('default'); % Initialize the random number generator stream
NeNodeB = 4;    % Number of eNodeB

% Create eNodeB configurations
enb = cell(1,NeNodeB);
for i=1:NeNodeB
   enb{i}=lteRMCDL('R.5');   % Get configuration based on RMC
   enb{i}.NCellID = i-1;     % Set unique cell identity
   enb{i}.TotSubframes = 1;  % Number of subframes to generate
   enb{i}.NPRSRB = 2;        % Bandwidth of PRS  in resource blocks 
   enb{i}.IPRS = 0;          % PRS configuration index
   enb{i}.PRSPeriod = 'On';  % PRS present in all subframes
   enb{i}.Position = hPositioningPosition(i-1, NeNodeB); % eNodeB position
end        

% Use the first eNodeB configuration for general settings
info = lteOFDMInfo(enb{1});

%% Plot Location of eNodeBs and UE    
% The positions of the eNodeBs and the UE are plotted for reference. The UE
% lies at (0,0) and the eNodeBs are distributed around the UE.

hPositioningPlotPositions(enb);

%% Generate Transmissions and Plot Transmitted Waveforms    
% For each eNodeB, a transmission is made consisting of solely the PRS. An
% empty resource grid is created and a PRS is generated and mapped onto the
% grid using <matlab:doc('ltePRS') ltePRS> and <matlab:doc('ltePRSIndices')
% ltePRSIndices>. The resultant grid is OFDM modulated to produce a
% transmit waveform.

tx = cell(1,NeNodeB);
for i = 1:NeNodeB
    grid = lteDLResourceGrid(enb{i});             % Empty resource grid
    grid(ltePRSIndices(enb{i})) = ltePRS(enb{i}); % Map PRS REs
    tx{i} = lteOFDMModulate(enb{i}, grid);        % OFDM modulate
end

% Plot transmit waveforms from each eNodeB
hPositioningPlotTx(tx, info.SamplingRate);

%% Compute delays from eNodeBs to UEs    
% Using the known eNodeB positions, the time delay from each eNodeB to the
% UE is calculated using the distance between the UE and eNodeB, |radius|,
% and the speed of propagation (speed of light). Using knowledge of the
% sampling rate, |info.SamplingRate|, the sample delay and stored in
% |sampleDelay|. These variables will be used to model the environment
% between the eNodeBs and the UE but the information will NOT be provided
% to the UE.

speedOfLight = 299792458.0; % Speed of light in m/s

sampleDelay = zeros(1, NeNodeB);
radius = cell(1, NeNodeB);
for i = 1:NeNodeB
   [~, radius{i}] = cart2pol(enb{i}.Position(1), enb{i}.Position(2));
   delay = radius{i}/speedOfLight;                  % Delay in seconds
   sampleDelay(i) = round(delay*info.SamplingRate); % Delay in samples 
end

%% Create Sum of Received Waveforms and Plot Received Waveform
% The received signal at the UE is modeled by delaying each eNodeB
% transmission according to the values in |sampleDelay|, and attenuating
% the received signal from each eNodeB using the values in |radius| in
% conjunction with an implementation of the TR36.814 [ <#10 1> ] Urban
% Macro Line Of Sight (LOS) path loss model. The received waveform from
% each eNodeB is padded with zeros to ensure all waveforms are the same
% length.

sumrx = zeros(length(tx{1})+max(sampleDelay), 1);
rx = cell(1, NeNodeB);
for i = 1:NeNodeB
    % Urban Macro LOS path loss per TR36.814
    PLdB = hPositioningPathLoss(radius{i}, 2.1e9);
    PL = 10^(PLdB/10);
    
    % Add delay, pad and attenuate
    rx{i} = [zeros(sampleDelay(i), 1); tx{i}; ... 
                zeros(max(sampleDelay)-sampleDelay(i), 1)]/ sqrt(PL);    
            
    % Sum waveforms from all eNodeBs
    sumrx = sumrx + rx{i};
end

% Plot received waveforms
hPositioningPlotRx(rx, info.SamplingRate);  

%% Estimate Arrival Times
% The arrival times of the signals from each eNodeB are established at the
% UE by correlating the incoming signal with a local PRS generated with the
% cell identity of each eNodeB. Note that the absolute arrival times cannot
% be used at the UE to calculate it's position as it has no knowledge of
% how far away the eNodeBs are, only the difference in distances given by
% the difference in arrival times. Therefore the peak correlation for each
% eNodeB is used as a delay estimate to allow comparison.

ref = cell(1, NeNodeB);
corr = cell(1, NeNodeB);
delayEst = zeros(1, NeNodeB);
for i = 1:NeNodeB        
    % Generate reference PRS
    grid = lteDLResourceGrid(enb{i});
    grid(ltePRSIndices(enb{i})) = ltePRS(enb{i});
    ref{i} = lteOFDMModulate(enb{i},grid);
    
    % Correlate received signal with each reference PRS
    c = abs(xcorr(sumrx,ref{i}));
    
    % Reduced length of correlation vector for positioning and plotting
    c(1:length(sumrx)) = [];    % Remove meaningless result at beginning
    corr{i} = c(1:info.Nfft);   % Extract an OFDM symbol's worth of data
    
    % Delay estimate is at point of maximum correlation
    delayEst(i) = find(corr{i}==max(corr{i}));
end

% Plot correlation
hPositioningPlotCorr(corr, info.SamplingRate);

%% Compute TDOA and plot constant TDOA hyperbolas
% Using the arrival times, the time differences of arrival between each
% pair of eNodeBs is calculated using <matlab:edit('hPositioningTDOA.m')
% hPositioningTDOA.m>. The time differences are used to compute hyperbolas
% of constant delay difference, which are plotted relative to the known
% eNodeB positions and intersect at the position of the UE.

% Estimate time difference of arrival from each eNodeB
tdoa = hPositioningTDOA(delayEst,info.SamplingRate);

% Plot hyperbolas
figure(1);
for j = 1:NeNodeB
    for i = (j+1):NeNodeB
        dd = tdoa(i,j)*speedOfLight; % Delay distance
        [x, y] = hPositioningTDOACurve(enb{i}.Position, ...
            enb{j}.Position, dd);
        plot(x, y, 'k:', 'LineWidth', 2);            
    end
end

% Plot eNodeB marker
plot(0, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'c', 'LineWidth', 2);

%% Appendix
% This example uses the helper functions:
%
% * <matlab:edit('hPositioningPosition.m') hPositioningPosition.m>
% * <matlab:edit('hPositioningPlotPositions.m')
% hPositioningPlotPositions.m>
% * <matlab:edit('hPositioningPlotTx.m') hPositioningPlotTx.m>
% * <matlab:edit('hPositioningPathLoss.m') hPositioningPathLoss.m>
% * <matlab:edit('hPositioningPlotRx.m') hPositioningPlotRx.m>
% * <matlab:edit('hPositioningPlotCorr.m') hPositioningPlotCorr.m>
% * <matlab:edit('hPositioningTDOA.m') hPositioningTDOA.m>
% * <matlab:edit('hPositioningTDOACurve.m') hPositioningTDOACurve.m>

%% Selected Bibliography
% # 3GPP TR36.814.

displayEndOfDemoMessage(mfilename)


