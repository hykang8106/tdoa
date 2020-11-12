function [] = tdoa_shell_2d(snr_db, sensor_number, Trials)
% where is crow?
% 2d version of TDOAshell.m
% original(TDOAshell.m): https://github.com/StevenJL/tdoa_localization
% sensor number must be greater than 3
%
% [input]
% - snr_db: snr in db
% - sensor_number: sensor number, must be greater than 3
% - Trials: trial number
% [usage]
% tdoa_shell_2d(0, 8, 10)

% #############################################
% #### N must be greater than 3, N >= 4
% #############################################
if sensor_number < 4
    fprintf('sensor number MUST be greater than 3\n');
    return;
end

% radius ratio between target and sensor. radius_ratio >= 1
radius_ratio = 1.5; % if radius_ratio > 1, source may exist outside sensor array

% redundant for readability
N = sensor_number; % sensor(microphone) number in array

sound_speed = 340.29;

[wave, fs] = audioread('crow_cry.wav');
fs;
wave = wave(:, 1);
scale = 0.8 / max(wave);
wave = scale * wave;

Radius = 100; % sensor(microphone) array radius

% generate sensor(microphone) array position(uniform circular array)
Theta = linspace(0, 2 * pi, N + 1);
X = Radius * cos(Theta(1 : end - 1));
Y = Radius * sin(Theta(1 : end - 1));
Sen_position = [X.', Y.']; % microphone array position

True_position = zeros(Trials, 2);
Est_position = zeros(Trials, 2);

% Generate position of source
for i = 1 : Trials
    r = rand * Radius * radius_ratio;
    t = rand * 2 * pi;
    True_position(i, 1) = r * cos(t);
    True_position(i, 2) = r * sin(t);
end
True_position;

% Generate distances between source and sensor
Distances = zeros(Trials, N);
for i = 1 : Trials
    for j = 1 : N
        x1 = True_position(i, 1);
        y1 = True_position(i, 2);
        
        x2 = Sen_position(j, 1);
        y2 = Sen_position(j, 2);
        
        Distances(i, j) = sqrt((x1 - x2)^2 + (y1 - y2)^2);
    end
end
Distances;
TimeDelay = Distances / sound_speed; % sound speed = 340.29 m/s
% what is padding?!
Padding = TimeDelay * fs;

% Generate the signals
mic = cell(1, N);
pre_length = zeros(1, N);
for i = 1 : Trials
    for j = 1 : N
        % prepend zero padding before wave sample
        % this make signal delayed when sensor receive signal
        mic{j} = [zeros(round(Padding(i, j)), 1); wave];
        pre_length(j) = length(mic{j});
    end
    
    m = max(pre_length);    
    c = m - pre_length;
    
    for j = 1 : N
        % append zero padding after wave sample
        mic{j} = [mic{j}; zeros(c(j), 1)];
        % attenuate signal
        mic{j} = mic{j} / (Distances(i, j)^2);
        % add awgn
        mic{j} = awgn(mic{j}, snr_db, 'measured', 'db'); 
    end
    
    multitrack = cell2mat(mic);
    size(multitrack)
    
    [x, y] = Locate(Sen_position, multitrack, fs, sound_speed);
    Est_position(i, 1) = x;
    Est_position(i, 2) = y;
end

delta_position = Est_position - True_position;
position_error = (sqrt(delta_position(:, 1).^2 + delta_position(:, 2).^2))';
position_error

figure('Position', [554 558 1018 420]); 

subplot(1, 2, 1);
sample_length = length(wave);
x = (0 : sample_length - 1) / fs;
plot(x, wave);
xlim([x(1) x(end)]);
xlabel('time in sec');
ylabel('signal magnitude');
title('tx waveform');
grid on;

subplot(1, 2, 2);
sample_length = size(multitrack, 1);
x = (0 : sample_length - 1) / fs;
plot(x, multitrack);
xlim([x(1) x(end)]);
xlabel('time in sec');
ylabel('signal magnitude');
title(sprintf('rx waveform: snr = %d db, sensor number = %d', snr_db, sensor_number));
grid on;

figure;
plot(Sen_position(:, 1), Sen_position(:, 2), 'go', ...
    True_position(:, 1), True_position(:, 2), 'bd', ...
    Est_position(:, 1), Est_position(:, 2), 'r+', ...
    'LineWidth', 2);
legend('Sensor Position', 'True Position', 'Estimated Position');
xlabel('X coordinate of target');
ylabel('Y coordinate of target');
title(sprintf('TDOA Localization: sensor = %d, snr = %d db, trial = %d', ...
    sensor_number, snr_db, Trials));
axis equal;
axis([-1 1 -1 1] * Radius * radius_ratio * 1.1);
grid on;

end

%%
function [x, y] = Locate(Sen_position, multitrack, fs, sound_speed)
% for tdoa linear equation, see "TDOA_Localization.pdf"(https://github.com/StevenJL/tdoa_localization)

len = size(Sen_position, 1);
timedelayvec = zeros(len, 1);
for i = 1 : len
    timedelayvec(i) = timedelayfunc(multitrack(:, 1), multitrack(:, i), fs);
end
timedelayvec

A = zeros(len, 1);
B = zeros(len, 1);
C = zeros(len, 1);

x1 = Sen_position(1, 1);
y1 = Sen_position(1, 2);

x2 = Sen_position(2, 1);
y2 = Sen_position(2, 2);
    
for i = 3 : len
    xi = Sen_position(i, 1);
    yi = Sen_position(i, 2);

    A(i) = (1/(sound_speed * timedelayvec(i))) * (-2 * x1 + 2 * xi) - ...
        (1/(sound_speed * timedelayvec(2))) * (-2 * x1 + 2 * x2);
    B(i) = (1/(sound_speed * timedelayvec(i))) * (-2 * y1 + 2 * yi) - ...
        (1/(sound_speed * timedelayvec(2))) * (-2 * y1 + 2 * y2);

    Sum1 = (x1^2) + (y1^2) - (xi^2) - (yi^2);
    Sum2 = (x1^2) + (y1^2) - (x2^2) - (y2^2);

    C(i) = sound_speed * (timedelayvec(i) - timedelayvec(2)) + ...
        (1/(sound_speed * timedelayvec(i))) * Sum1 - ...
        (1/(sound_speed * timedelayvec(2))) * Sum2;
end
A;
B;
C;

M = [A(3 : end) B(3 : end)];
C = C(3 : end);

Minv = pinv(M);
T = Minv * (-C);
x = T(1);
y = T(2);

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

