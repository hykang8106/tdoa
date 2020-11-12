function [] = tdoa_shell_2d(snr_db)
% 2d version of TDOAshell.m
% original(TDOAshell.m): https://github.com/StevenJL/tdoa_localization
% sensor number must be greater than 3
%
% [input]
% -snr_db: snr in db
% [usage]
% tdoa_shell_2d(0)

sound_speed = 340.29;

[wave, fs] = audioread('crow_cry.wav');
fs;
wave = wave(:, 1);
scale = 0.8 / max(wave);
wave = scale * wave;

Trials = 10;

Radius = 50; % sensor(microphone) array radius
% #############################################
% #### N must be greater than 3, N >= 4
% #############################################
N = 8; % sensor(microphone) number in array

% generate sensor(microphone) array position(uniform circular array)
Theta = linspace(0, 2 * pi, N + 1);
X = Radius * cos(Theta(1 : end - 1));
Y = Radius * sin(Theta(1 : end - 1));
Sen_position = [X.', Y.']; % microphone array position

True_position = zeros(Trials, 2);
Est_position = zeros(Trials, 2);

% Generate position of source
for i = 1 : Trials
    r = rand * Radius;
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
for i = 1 : Trials
    % prepend zero padding before wave sample
    % this make signal delayed when sensor receive signal
    mic1 = [zeros(round(Padding(i, 1)), 1); wave];
    mic2 = [zeros(round(Padding(i, 2)), 1); wave];
    mic3 = [zeros(round(Padding(i, 3)), 1); wave];
    mic4 = [zeros(round(Padding(i, 4)), 1); wave];
    mic5 = [zeros(round(Padding(i, 5)), 1); wave];
    mic6 = [zeros(round(Padding(i, 6)), 1); wave];
    mic7 = [zeros(round(Padding(i, 7)), 1); wave];
    mic8 = [zeros(round(Padding(i, 8)), 1); wave];
    
    l1 = length(mic1);
    l2 = length(mic2);
    l3 = length(mic3);
    l4 = length(mic4);
    l5 = length(mic5);
    l6 = length(mic6);
    l7 = length(mic7);
    l8 = length(mic8);
    
    m = max([l1 l2 l3 l4 l5 l6 l7 l8]);
    
    c = [m - l1, m - l2, m - l3, m - l4, m - l5, m - l6, m - l7, m - l8];
    
    % append zero padding after wave sample
    mic1 = [mic1; zeros(c(1), 1)];
    mic2 = [mic2; zeros(c(2), 1)];
    mic3 = [mic3; zeros(c(3), 1)];
    mic4 = [mic4; zeros(c(4), 1)];
    mic5 = [mic5; zeros(c(5), 1)];
    mic6 = [mic6; zeros(c(6), 1)];
    mic7 = [mic7; zeros(c(7), 1)];
    mic8 = [mic8; zeros(c(8), 1)];
    
    % attenuate signal and add awgn
    mic1 = mic1 / Distances(i, 1);  mic1 = awgn(mic1, snr_db, 'measured', 'db'); 
    mic2 = mic2 / Distances(i, 2);  mic2 = awgn(mic2, snr_db, 'measured', 'db'); 
    mic3 = mic3 / Distances(i, 3);  mic3 = awgn(mic3, snr_db, 'measured', 'db'); 
    mic4 = mic4 / Distances(i, 4);  mic4 = awgn(mic4, snr_db, 'measured', 'db'); 
    mic5 = mic5 / Distances(i, 5);  mic5 = awgn(mic5, snr_db, 'measured', 'db'); 
    mic6 = mic6 / Distances(i, 6);  mic6 = awgn(mic6, snr_db, 'measured', 'db'); 
    mic7 = mic7 / Distances(i, 7);  mic7 = awgn(mic7, snr_db, 'measured', 'db'); 
    mic8 = mic8 / Distances(i, 8);  mic8 = awgn(mic8, snr_db, 'measured', 'db'); 
    
    multitrack = [mic1, mic2, mic3, mic4, mic5, mic6, mic7, mic8];
    size(multitrack);
    
    [x, y] = Locate(Sen_position, multitrack, fs, sound_speed);
    Est_position(i, 1) = x;
    Est_position(i, 2) = y;
end

delta_position = Est_position - True_position;
position_error = (sqrt(delta_position(:, 1).^2 + delta_position(:, 2).^2))';
position_error

figure;
plot(multitrack);

figure;
plot(Sen_position(:, 1), Sen_position(:, 2), 'go', ...
    True_position(:, 1), True_position(:, 2), 'bd', ...
    Est_position(:, 1), Est_position(:, 2), 'r+', ...
    'LineWidth', 2);
legend('Sensor Position', 'True Position', 'Estimated Position');
xlabel('X coordinate of target');
ylabel('Y coordinate of target');
title('TDOA Hyperbolic Localization');
axis equal;
axis([-1 1 -1 1] * Radius * 1.1);
grid on;

end

%%
function [x, y] = Locate(Sen_position, multitrack, fs, sound_speed)

len = size(Sen_position, 1);
timedelayvec = zeros(len, 1);
for i = 1 : len
    timedelayvec(i) = timedelayfunc(multitrack(:, 1), multitrack(:, i), fs);
end
timedelayvec;

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

end

