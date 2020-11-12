function [] = ...
    batch_simulate_tdoa_sensor_fixed(sensor_length, sensor_exist_at_uca_center, uca_radius_meter, ...
    snr_db, radius_ratio, sensor_position_filename, trial_length, delta_distance, target_signal_spec)
% batch version of simulate_tdoa_sensor_fixed.m
% i want to see how position error is varied according to target position and sensor length/shape
%
% [usage]
% batch_simulate_tdoa_sensor_fixed(5, 1, 5e3, 5, 1.5, [], 100, [100, 100], [15,2,1])
% batch_simulate_tdoa_sensor_fixed(5, 1, 5e3, 5, 1, 'sensor_position.txt', 300, [100, 100], [15,2,1])

% if sensor_length > 7 || sensor_length < 3
%     fprintf('sensor number: 3 ~ 7\n');
%     return;
% end

% ### MUST USE [], not implemented yet
% sensor_position_filename = [];
% sensor_position_filename = '';

% ############################
% #### MUST be empty!!!
% ############################
% fading cant be applied in batch_simulate_tdoa_sensor_fixed
% because batch_simulate_tdoa_sensor_fixed show position error versus all target position
% target position near sensor is not good to apply fading
rician_param = [];

plot_signal = 0;
plot_position = 0;
use_only_torrieri_method = 1;

fs = 3.84e6;
overlap_criterion_meter = 1;
% overlap_criterion_meter = 10;
% overlap_criterion_meter = physconst('LightSpeed') / fs / 2;

delta_x = delta_distance(1);
delta_y = delta_distance(2);

% delta_x = 100;
% delta_y = 100;

% for saving result
filename = 'tdoa_result_sensor_fixed.mat';

if ~isempty(sensor_position_filename)
    [sensor_position, uca_radius_meter] = get_sensor_position_from_file(sensor_position_filename);
    sensor_length = size(sensor_position, 1);
%     sensor_position = load(sensor_position_filename);
else
    sensor_position = get_sensor_position(sensor_length, uca_radius_meter, sensor_exist_at_uca_center);
end

if sensor_length > 7 || sensor_length < 3
    fprintf('sensor number: 3 ~ 7\n');
    return;
end

target_radius_meter = uca_radius_meter * radius_ratio;

x_range = [-1, 1] * target_radius_meter;
y_range = [-1, 1] * target_radius_meter;

x = x_range(1) : delta_x : x_range(end);
x_length = length(x);
y = y_range(1) : delta_y : y_range(end);
y_length = length(y);

error_torrieri = nan(x_length, y_length);
cep_mean = nan(x_length, y_length);
gdop_mean = nan(x_length, y_length);

h = waitbar(0, 'Please wait...', 'name', 'batch program progress');
steps = x_length * y_length;

u = tic;

for i = 1 : x_length
    for j = 1 : y_length
        
        position_error_torrieri = zeros(1, trial_length);
        % ###### for computing position estimation covariance
        x_est_torrieri = zeros(2, trial_length);
        cep = zeros(1, trial_length);
        gdop = zeros(1, trial_length);
        
        target_position = [x(i), y(j)]
        
        if ~check_target_position_is_good(target_position, target_radius_meter)
            fprintf('###### target position is out of range\n');
            return;
        end
        
        if check_target_overlap_sensor(sensor_position, target_position, overlap_criterion_meter)
            continue; % target overlap sensor, so skip
        end
        
        for k = 1 : trial_length
            [position_error_torrieri(k), x_est_torrieri(:, k), cep(k), gdop(k), ...
                lambda1, lambda2, theta_rad] = ...
                sub_simulate_tdoa_sensor_fixed(sensor_position, snr_db, target_position, ...
                plot_signal, plot_position, use_only_torrieri_method, target_radius_meter, rician_param, ...
                target_signal_spec);
        end
        
        position_error_torrieri;
        mean_torrieri = mean(position_error_torrieri);
        std_torrieri = std(position_error_torrieri);
        
        error_torrieri(i, j) = mean_torrieri;
        cep_mean(i, j) = mean(cep);
        gdop_mean(i, j) = mean(gdop);
        
        step = y_length * (i - 1) + j;
        w = toc(u);
        % reference: 
        % http://stackoverflow.com/questions/12210583/is-there-a-matlab-function-to-convert-elapsed-seconds-to-hhmmss-format
        z = fix(mod(w, [0, 3600, 60])./[3600, 60, 1]);
        waitbar(step / steps, h, sprintf('%f %%, elapsed time = %d : %02d : %02d', ...
            step / steps * 100, z(1), z(2), z(3)));
        
    end
end

close(h);

error_torrieri
cep_mean
gdop_mean

% #### append date string to filename
[p, name, e] = fileparts(filename);
filename = [name, '_', datestr(now, 'yymmddHHMM'), e];
save(filename, 'sensor_position', 'snr_db', 'trial_length', 'x', 'y', 'target_radius_meter', ...
    'delta_x', 'delta_y', 'error_torrieri', 'cep_mean', 'gdop_mean');

% datestr(now, 'yymmddHHMMSS')

% position_error_torrieri;
% mean_torrieri = mean(position_error_torrieri)
% std_torrieri = std(position_error_torrieri)

% sensor_vec_len = length(sensor_length);
% snr_vec_len = length(snr_db);
% 
% position_error_torrieri = zeros(trial_length, snr_vec_len, sensor_vec_len);
% position_error_linear = zeros(trial_length, snr_vec_len, sensor_vec_len);
% position_error_hyperbolic = zeros(trial_length, snr_vec_len, sensor_vec_len);
% 
% for i = 1 : sensor_vec_len
%     for j = 1 : snr_vec_len
%         for k = 1 : trial_length
%             [position_error_torrieri] = ...
%                 simulate_tdoa_sensor_fixed(sensor_length, uca_radius_meter, snr_db, target_position, ...
%                 sensor_exist_at_uca_center, radius_ratio)
%             
%             [error_hyperbolic, ...
%                 error_linear, ...
%                 position_error_torrieri(k, j, i)] = ...
%                 simulate_tdoa_target_fixed(sensor_length(i), snr_db(j), randomize_sensor_distance);
%             
% %             if ~isempty(error_hyperbolic)
% %                 position_error_hyperbolic(k, j, i) = error_hyperbolic;
% %             end
% %             
% %             if ~isempty(error_linear)
% %                 position_error_linear(k, j, i) = error_linear;
% %             end
%         end
%     end
% end
% 
% position_error_torrieri;
% position_error_hyperbolic;
% 
% torrieri_mean = squeeze(mean(position_error_torrieri))
% torrieri_std = squeeze(std(position_error_torrieri))
% 
% % hyperbolic is valid when sensor length is only 3
% idx = find(sensor_length == 3, 1);
% if ~isempty(idx)
%     hyperbolic_mean = mean(position_error_hyperbolic(:, :, idx))
%     hyperbolic_std = std(position_error_hyperbolic(:, :, idx))
% end
% 
% % linear is valid when sensor length is greater than 3
% idx = find(sensor_length >= 4);
% if ~isempty(idx)
%     linear_sensor_len = sensor_length(idx);
%     
%     position_error_linear = position_error_linear(:, :, idx);
%     position_error_linear;
%     [trial_len, snr_len, sensor_len] = size(position_error_linear);
%     
%     % ###### special care is needed because linear have nan ########
%     linear_mean = zeros(snr_len, sensor_len);
%     linear_std = zeros(snr_len, sensor_len);
%     
%     for n = 1 : sensor_len
%         for m = 1 : snr_len
%             tmp = squeeze(position_error_linear(:, m, n));
%             idx = ~isnan(tmp);
%             tmp = tmp(idx);
%             
%             linear_mean(m, n) = mean(tmp);
%             linear_std(m, n) = std(tmp);
%         end
%     end  
%     
%     linear_mean
%     linear_std
% end
% 
% save(filename, 'sensor_length', 'snr_db', 'trial_length', ...
%     'torrieri_mean', 'torrieri_std', ...
%     'hyperbolic_mean', 'hyperbolic_std', ...
%     'linear_mean', 'linear_std');

end

