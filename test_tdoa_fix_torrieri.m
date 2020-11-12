function [] = test_tdoa_fix_torrieri(rf_sensor_length, snr_db, randomize_sensor_distance)
% #### test tdoa_fix_torrieri.m, copied from tdoa_position.m ####
%
% target location using tdoa method
% modified to simulate tdoa location of target emitter used for spectrum monitoring
% original was PRSPositioningExample.m in LTE System Toolbox example
% target emitter is located at (0,0)
% 
% [input]
% - rf_sensor_length: number of rf sensor, 3 ~ 5
%   ofcom experience show sensor number greater than 5 dont improve location accuracy
% - snr_db: snr in db
% - randomize_sensor_distance: 0 or 1, 
%   ############# temporary set = 0(reference sensor is sensor 1)
% [usage]
% test_tdoa_fix_torrieri(4, 0, 0)
%

%%
if rf_sensor_length > 7 || rf_sensor_length < 3
    fprintf('rf sensor number: 3 ~ 7\n');
    return;
end
rf_sensor = cell(1,rf_sensor_length);

plot_signal = 0; % control plot of tx, rx, correlation
plot_position = 0; % control plot of target, sensor, tdoa curve

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

% tx = cell(1,1);
grid = lteDLResourceGrid(enb);          % Empty resource grid
grid(ltePRSIndices(enb)) = ltePRS(enb); % Map PRS REs
tx{1} = lteOFDMModulate(enb, grid);        % OFDM modulate

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
if plot_signal
    my_hPositioningPlotCorr(corr, info.SamplingRate, ref_sensor);
end

%% Compute TDOA between reference sensor and other sensor, and plot constant TDOA hyperbolas

delayEst
delta_distance_meter = delayEst / info.SamplingRate * speedOfLight

if plot_position
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
    plot(0, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'r', 'LineWidth', 2);
    % plot(0, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'c', 'LineWidth', 2);
end

%%
locate_target_tdoa(rf_sensor, ref_sensor, delta_distance_meter, rx, info.SamplingRate);

end

%%
function [] = locate_target_tdoa(rf_sensor, ref_sensor_number, delta_distance_meter, rx, fs)
% locate target using tdoa method
%
% ### location method ###
% (1) hyperbolic equation solver
% (2) linear equation solver
% (3) torrieri algorithm
%
% [input]
% - rf_sensor: cell array having rf sensor position
% - ref_sensor: reference sensor number
% - delta_distance_meter: 1 x sensor_length
% - rx: cell array having rx signal, only for linear_eq_solver function
%   ##### rx input is NOT nice, but this is temporary for linear equation solver
%   ##### rx input is used for computing time delay using cross correlation of rx signal
% - fs: sample rate

rf_sensor;
sensor_length = length(rf_sensor);

% convert cell to sensor position matrix(dimension = sensor_length x 2)
sensor_position = zeros(sensor_length, 2);
for n = 1 : sensor_length
    rf_sensor{n}.Position;
    sensor_position(n, :) = rf_sensor{n}.Position;
end
sensor_position

if sensor_length == 3
    % hyperbolic equation solver
    sensor_1 = sensor_position(1, :);
    sensor_2 = sensor_position(2, :);
    sensor_3 = sensor_position(3, :);
    % initial guess: mean of sensor position
    x0 = mean(sensor_position);
    
    delta_distance_12 = delta_distance_meter(2);
    delta_distance_13 = delta_distance_meter(3);
    [x_est, fval, exitflag] = ...
        hyperbolic_eq_solver(x0, sensor_1, sensor_2, sensor_3, delta_distance_12, delta_distance_13);
    if exitflag == 1
        x_est
        % #################################
        % ### exact solution is [0 0]
        % #################################
        position_error = sqrt(x_est(1)^2 + x_est(2)^2);
        fprintf('##### hyperbolic equation solver: position error = %f meter\n', position_error);
    else
        fprintf('############### fsolve error(NOT converged to a root): exit flag = %d\n', exitflag);
    end
    % eval('[x_ref, fval, exitflag] = hyperbolic_eq_solver(x0, sensor_1, sensor_2, sensor_3, delta_distance_12, delta_distance_13)');
    
    % torrieri algorithm
    x_est = thanks_torrieri(sensor_position, ref_sensor_number, delta_distance_meter);
    position_error = sqrt(x_est(1)^2 + x_est(2)^2);
    fprintf('##### torrieri algorithm: position error = %f meter\n', position_error);
elseif sensor_length >= 4
    % linear equation solver
    x_est = linear_eq_solver(sensor_position, ref_sensor_number, rx, fs)
    
    % #################################
    % ### exact solution is [0 0]
    % #################################
    position_error = sqrt(x_est(1)^2 + x_est(2)^2);
    fprintf('##### linear equation solver: position error = %f meter\n', position_error);
    
    % torrieri algorithm
    x_est = thanks_torrieri(sensor_position, ref_sensor_number, delta_distance_meter);
    position_error = sqrt(x_est(1)^2 + x_est(2)^2);
    fprintf('##### torrieri algorithm: position error = %f meter\n', position_error);
else
    % ### NOT reach here ###
    error('number of sensors must be greater than 2\n');
end

end

%%
function [x, fval, exitflag] = ...
    hyperbolic_eq_solver(x0, sensor_1, sensor_2, sensor_3, delta_distance_12, delta_distance_13)

options = optimoptions('fsolve','Display','off'); % Option to display output
f = @(x)hyperbolic_fun(x, sensor_1, sensor_2, sensor_3, delta_distance_12, delta_distance_13);
[x, fval, exitflag] = fsolve(f,x0,options); % Call solver

end

%%
function F = hyperbolic_fun(x, sensor_1, sensor_2, sensor_3, delta_distance_12, delta_distance_13)

F = [sqrt((sensor_1(1) - x(1))^2 + (sensor_1(2) - x(2))^2) - ...
    sqrt((sensor_2(1) - x(1))^2 + (sensor_2(2) - x(2))^2) + ...
    delta_distance_12;
    sqrt((sensor_1(1) - x(1))^2 + (sensor_1(2) - x(2))^2) - ...
    sqrt((sensor_3(1) - x(1))^2 + (sensor_3(2) - x(2))^2) + ...
    delta_distance_13];

end

%%
function [T] = linear_eq_solver(sensor_position, ref_sensor, rx, fs)
% for tdoa linear equation, see "TDOA_Localization.pdf"(https://github.com/StevenJL/tdoa_localization)
%
% [input]
% - sensor_position: sensor position, dimension = sensor_length x 2
% - ref_sensor: reference sensor number
% - rx: cell array having rx signal

v = physconst('LightSpeed');

% rx was cell array, so convert matrix
rx = cell2mat(rx);
size(rx);

len = size(sensor_position, 1);
timedelayvec = zeros(len, 1);
for i = 1 : len
    timedelayvec(i) = timedelayfunc(rx(:, 1), rx(:, i), fs);
end
timedelayvec;

A = zeros(len, 1);
B = zeros(len, 1);
C = zeros(len, 1);

x1 = sensor_position(1, 1);
y1 = sensor_position(1, 2);

x2 = sensor_position(2, 1);
y2 = sensor_position(2, 2);
    
for i = 3 : len
    xi = sensor_position(i, 1);
    yi = sensor_position(i, 2);

    A(i) = (1/(v * timedelayvec(i))) * (-2 * x1 + 2 * xi) - ...
        (1/(v * timedelayvec(2))) * (-2 * x1 + 2 * x2);
    B(i) = (1/(v * timedelayvec(i))) * (-2 * y1 + 2 * yi) - ...
        (1/(v * timedelayvec(2))) * (-2 * y1 + 2 * y2);

    Sum1 = (x1^2) + (y1^2) - (xi^2) - (yi^2);
    Sum2 = (x1^2) + (y1^2) - (x2^2) - (y2^2);

    C(i) = v * (timedelayvec(i) - timedelayvec(2)) + ...
        (1/(v * timedelayvec(i))) * Sum1 - ...
        (1/(v * timedelayvec(2))) * Sum2;
end
A;
B;
C;

M = [A(3 : end) B(3 : end)];
C = C(3 : end);

Minv = pinv(M);
T = Minv * (-C);
% x = T(1);
% y = T(2);

end

%%
function out = timedelayfunc(x, y, fs)

c = xcorr(x, y);
[C, I] = max(c);
out = ((length(c) + 1)/2 - I)/fs;
% (length(c) + 1)/2 - I

% I - length(x)
% c(1:length(x)) = [];
% 
%     % Reduced length of correlation vector for positioning and plotting
%     c(1:length(ref_signal)) = [];    % Remove meaningless result at beginning
%     corr{i} = c(1:info.Nfft);   % Extract an OFDM symbol's worth of data
%     
%     % Delay estimate is at point of maximum correlation
%     delayEst(i) = find(corr{i} == max(corr{i}));

end

%%
function [x_est] = thanks_torrieri(sensor_position, ref_sensor_number, delta_distance_meter)
%
% [input]
% - sensor_position: sensor_length x 2
% - ref_sensor: reference sensor number
% - delta_distance_meter: 1 x sensor_length
% [output]
% - x_est: dimension = 2 x 1

iteration_length = 5;
% criterion for iteration stop: 
% when absolute difference between initial and estimated is greater than criterion, iteration stop
iteration_criterion = 10;

c = physconst('LightSpeed');

% % set tdoa_position.m input:
% % randomize_sensor_distance = 0 (sensor 1 is always reference sensor, emitter is nearest to sensor 1)
% % rf_sensor_length = 4
% sensor_1 = [4.8541   -0.7114]*1e3;
% sensor_2 = [3.9030    4.7709]*1e3;
% sensor_3 = [-6.3239    1.8785]*1e3;
% sensor_4 = [-4.4540   -7.0001]*1e3;
% 
% delta_distance_12 = 1.2491e3;
% delta_distance_13 = 1.7176e3;
% delta_distance_14 = 3.3571e3;

% x0 = [0; 0];  % Make a starting guess at the solution
% 
% [x_ref, fval, exitflag] = hyperbolic_eq_solver(x0, sensor_1, sensor_2, sensor_3, delta_distance_12, delta_distance_13)
% % eval('[x_ref, fval, exitflag] = hyperbolic_eq_solver(x0, sensor_1, sensor_2, sensor_3, delta_distance_12, delta_distance_13)');
% 
% % ######## considering average of hyperbolic equation solver output #######
% [x_ref, fval, exitflag] = hyperbolic_eq_solver(x0, sensor_1, sensor_3, sensor_4, delta_distance_13, delta_distance_14)

% #### consider average of sensor position
% ############# this option is most reasonable ###########
% x_ref = mean([sensor_1; sensor_2; sensor_3; sensor_4])'
x_ref = mean(sensor_position)';

% ##### need minus sign #####
% Ht = -[delta_distance_12, delta_distance_13, delta_distance_14]' / c;
Ht = -delta_distance_meter(2 : end)' / c; % dont forget to remove first element of delta_distance_meter

% % sensor dimension: 2 x sensor_length
% sensor = [sensor_1' sensor_2' sensor_3' sensor_4'];

% x_ref: column vector (2 x 1)
Rhat = zeros(2, iteration_length);
for n = 1 : iteration_length
    [x_est, lambda1, lambda2, cep, gdop, theta_rad] = torrieri(x_ref, sensor_position', Ht);
    Rhat(:, n) = x_est;
    
    % test iteration criterion
    if abs(sqrt(x_est(1)^2 + x_est(2)^2) - sqrt(x_ref(1)^2 + x_ref(2)^2)) <= iteration_criterion
        break;
    end
    
    x_ref = x_est;
end
Rhat
x_est, lambda1, lambda2, cep, gdop, theta_rad
position_error = sqrt(x_est(1)^2 + x_est(2)^2);

end

%%
function [Rhat, lambda1, lambda2, cep, gdop, theta_rad] = torrieri(R0, sensor, Ht)
% thanks torrieri
% see "Statistical Theory of Passive Location Systems" paper written by torrieri
%
% [input]
% - R0: initial target position guess, 2 x 1
% - sensor: sensor position, 2 x sensor_length
% - Ht: tdoa between reference sensor and other sensor, (sensor_length - 1) x 1
%   Ht means H * t in eq (81) where t is left term in eq (79) and H is eq (82)
%   in current version, reference sensor is always sensor 1 
% [output]
% - Rhat: estimated target position, 2 x 1
% - lambda1: major axis length of concentration ellipse
% - lambda2: minor axis length of concentration ellipse
% - cep: circular error probable
% - gdop: geometric dilution of precision
% - theta_rad: coordinate rotation angle in rad

c = physconst('LightSpeed');
fs = 3.84e6;
sp = 1 / fs; % sample period
% sp = 1e-7 * sp

bw = 1e6; % signal bandwidth
snr_db = 10;

% see eq (94)
beta_r2 = pi^2 * bw^2 / 3;

% ##### atv: arrival time variance, sigma_t2
% ##### how to compute E and N0? how to related with snr_db?
% see eq (91)
% atv = 1 / ((2 * E / N0) * beta_r2)
% E: energy in received signal, N0/2: two-sided noise power spectral denisty,
% beta_r2: function of bandwidth of signal, see eq (94) for communication signals
% N0 = kT where k is boltzmann constant, T is system noise temperature

% matlab command: systemp, noisepow

% k = physconst('Boltzmann')

% it is nonsense that we set arrival time variance to sample peariod. sample synch must be assumed
% atv = sp;

% (2 * E / N0) = 10
% atv = 1 / (10 * beta_r2);

% ### temporary setting: E/No is assumed to be equal to Eb/No
Eb_over_N0_db = 10;
Eb_over_N0 = 10^(Eb_over_N0_db / 10);
atv = 1 / (2 * Eb_over_N0 * beta_r2);

% get sensor length
sensor_length = size(sensor, 2);

% R0 = repmat(R0, 1, sensor_length);
% repeat R0 column vector to make matrix whose dimension is (2 x sensor_length)
F = (repmat(R0, 1, sensor_length) - sensor)'; % F dimension: sensor_length x 2

D0 = sqrt(F(:, 1).^2 + F(:, 2).^2);
% see eq (78)
F = F ./ repmat(D0, 1, 2);

% for n = 1 : sensor_length
%     F(n, :) = F(n, :) / sqrt(F(n, 1)^2 + F(n, 2)^2);
% end
% F

% ### valid only when sensor 1 is reference sensor
% see eq (82) (modified so that ref sensor is sensor 1)
% H = [...
%     1 -1 0 0;
%     1 0 -1 0;
%     1 0 0 -1;
%     ];
% % #### general form for arbitrary sensor length, but reference sensor MUST be sensor 1
% H = diag(-ones(1, sensor_length - 1), 1);
% H = H(1 : end - 1, :); % remove last row
% H(:, 1) = 1; % set first column to 1

% #### use make_torrieri_H_matrix.m #####
ref_sensor_number = 1;
[H] = make_torrieri_H_matrix(sensor_length, ref_sensor_number);

% atv(sigmat) is related to epsilon(eq (75), (76)), 
% which is arrival time measurement error accounts for
% propagation anomalies, receiver noise, and errors in the assumed station positions.
% see eq (96) (modified so that sensor number is 4)
% Ne = [...
%     atv 0 0 0;
%     0 atv 0 0;
%     0 0 atv 0;
%     0 0 0 atv;
%     ];
Ne = diag(ones(1, sensor_length) * atv); % general form for arbitrary sensor length

% see eq (85)
N = H * Ne * H';

FtHtNi = F' * H' * inv(N);

% see eq (86)
Rhat = R0 + c * inv(FtHtNi * H * F) * FtHtNi * (Ht - H * D0 / c);

% #### P: covariance matrix of estimated target position
% see eq (16),(87)
P = (c^2) * inv(FtHtNi * H * F);

E = eig(P);

% see eq (54)
lambda1 = .5 * [P(1,1) + P(2,2) + sqrt((P(1,1) - P(2,2))^2 + 4 * P(1,2)^2)];
% see eq (55)
lambda2 = .5 * [P(1,1) + P(2,2) - sqrt((P(1,1) - P(2,2))^2 + 4 * P(1,2)^2)];

% GDOP(geometric dilution of precision) is defined as 
% the ratio of the root-mean-square position error Er to the root mean-square ranging error.
% The GDOP indicates how much the fundamental ranging error is magnified 
% by the geometric relation among the transmitter position and the stations.
% see eq (89)
gdop = sqrt(trace(P)) / (c * sqrt(atv));

% CEP(circular error probable) is defined as the radius of the circle 
% that has its center at the mean and contains half the realizations of the random vector.
% CEP is a measure of the uncertainty in the location estimator x relative to its mean E[x].
% see eq (90)
cep = .75 * (c * sqrt(atv)) * gdop;
% see eq (74)
cep = .75 * sqrt(P(1,1) + P(2,2));

% coordinate rotation angle in radian, see eq (57)
theta_rad = .5 * atan(2 * P(1,2) / (P(1,1) - P(2,2)));
% coordinate rotation angle in degree, see eq (57)
theta_deg = .5 * atand(2 * P(1,2) / (P(1,1) - P(2,2)));

% concentration ellipse corresponding to probability Pe is defined 
% to be the particular ellipse for which Pe is the probability that x lies inside it. 
% Thus the concentration ellipse is a 2 dimensional measure of accuracy for an unbiased estimator.

% concentration ellipse corresponding to probability Pe: 
% major axes length = 2 * sqrt(k * lambda1) 
% minor axes length = 2 * sqrt(k * lambda2)
% where k = -2 * log(1 - Pe)
% see eq (60)

% for target position bias, use eq (84), (15)
% see eq (15)
% when linearization error is neglecting, bias is:
% b = inv(FtHtNi * H * F / c^2) * FtHtNi * E[n] 
% where E[n] is expectation of n which is arrival time measurement error
% see eq (75), (76) and see eq (107), (108)

% estimator target position(Rhat) accuracy:
% rms error, epsilon_r
% epsilon_r^2 = trace(P) + sum(bias.^2)
% where P is covariance matrix of estimated target position
% see eq (51), (52)

end



