function [filename] = batch_simulate_tdoa_target_fixed(sensor_length, snr_db, trial_length)
% batch version of simulate_tdoa_target_fixed.m
% ###################################################################################################
% #### how location error is varied according to sensor length, snr for each location method?
% #### location method: least square estimator, linear equation solver, hyperbolic equation solver
% ###################################################################################################
%
% [input]
% - sensor_length: sensor length. can be vector
% - snr_db: snr in db. can be vector
% - trial_length: trial length
%
% [usage]
% batch_simulate_tdoa_target_fixed([3:7], [-10:5:15], 100)

% ##############################################
% ### must set to 0: use all location method
% ##############################################
use_only_torrieri_method = 0;
randomize_sensor_distance = 1;

% ### must be zero. otherwise many figure is generated
plot_position = 0;
plot_signal = 0;

% for saving result
filename = 'tdoa_result_target_fixed.mat';

torrieri_mean = []; torrieri_std = [];
hyperbolic_mean = []; hyperbolic_std = [];
linear_mean = []; linear_std = [];

sensor_vec_len = length(sensor_length);
snr_vec_len = length(snr_db);

position_error_torrieri = zeros(trial_length, snr_vec_len, sensor_vec_len);
position_error_linear = zeros(trial_length, snr_vec_len, sensor_vec_len);
position_error_hyperbolic = zeros(trial_length, snr_vec_len, sensor_vec_len);

h = waitbar(0, 'Please wait...', 'name', 'batch program progress', 'windowstyle', 'modal');
steps = sensor_vec_len * snr_vec_len * trial_length;

u = tic;

for i = 1 : sensor_vec_len
    for j = 1 : snr_vec_len
        for k = 1 : trial_length
            [error_hyperbolic, ...
                error_linear, ...
                position_error_torrieri(k, j, i)] = ...
                simulate_tdoa_target_fixed(sensor_length(i), snr_db(j), randomize_sensor_distance, ...
                use_only_torrieri_method, plot_position, plot_signal);
            
            if ~isempty(error_hyperbolic)
                position_error_hyperbolic(k, j, i) = error_hyperbolic;
            end
            
            if ~isempty(error_linear)
                position_error_linear(k, j, i) = error_linear;
            end
            
            step = snr_vec_len * trial_length * (i - 1) + trial_length * (j - 1) + k;
            w = toc(u);
            % reference:
            % http://stackoverflow.com/questions/12210583/is-there-a-matlab-function-to-convert-elapsed-seconds-to-hhmmss-format
            z = fix(mod(w, [0, 3600, 60])./[3600, 60, 1]);
            waitbar(step / steps, h, sprintf('%.0f %%, elapsed time = %d : %02d : %02d', ...
                step / steps * 100, z(1), z(2), z(3)));
        end
    end
end

close(h);

position_error_torrieri;
position_error_hyperbolic;

torrieri_mean = squeeze(mean(position_error_torrieri))
torrieri_std = squeeze(std(position_error_torrieri))

% hyperbolic is valid when sensor length is only 3
idx = find(sensor_length == 3, 1);
if ~isempty(idx)
    hyperbolic_mean = mean(position_error_hyperbolic(:, :, idx))
    hyperbolic_std = std(position_error_hyperbolic(:, :, idx))
end

% linear is valid when sensor length is greater than 3
idx = find(sensor_length >= 4);
if ~isempty(idx)
    linear_sensor_len = sensor_length(idx);
    
    position_error_linear = position_error_linear(:, :, idx);
    position_error_linear;
    [trial_len, snr_len, sensor_len] = size(position_error_linear);
    
    % ###### special care is needed because linear equation solver have nan ########
    linear_mean = zeros(snr_len, sensor_len);
    linear_std = zeros(snr_len, sensor_len);
    
    for n = 1 : sensor_len
        for m = 1 : snr_len
            tmp = squeeze(position_error_linear(:, m, n));
            idx = ~isnan(tmp);
            tmp = tmp(idx);
            
            linear_mean(m, n) = mean(tmp);
            linear_std(m, n) = std(tmp);
        end
    end  
    
    linear_mean
    linear_std
end

[p, name, e] = fileparts(filename);
filename = [name, '_', datestr(now, 'yymmddHHMM'), e];
save(filename, 'sensor_length', 'snr_db', 'trial_length', ...
    'torrieri_mean', 'torrieri_std', ...
    'hyperbolic_mean', 'hyperbolic_std', ...
    'linear_mean', 'linear_std');

end




