function [] = ...
    batch_simulate_tdoa_sensor_fixed(sensor_length, sensor_exist_at_uca_center, uca_radius_meter, snr_db, ...
    radius_ratio, trial_length)
% batch version of simulate_tdoa_sensor_fixed.m
% i want to see how position error is varied according to target position and sensor length/shape
%
% [usage]
% batch_simulate_tdoa_sensor_fixed(5, 1, 6e3, 5, 1.0, 300)

% ### MUST 
USE [], not implemented yet
sensor_position_filename = [];
% sensor_position_filename = '';

fs = 3.84e6;
overlap_criterion_meter = physconst('LightSpeed') / fs; 

delta_x = 1000;
delta_y = 1000;

x_range = [-1, 1] * uca_radius_meter * radius_ratio;
y_range = [-1, 1] * uca_radius_meter * radius_ratio;

x = x_range(1) : delta_x : x_range(end);
x_length = length(x);
y = y_range(1) : delta_y : y_range(end);
y_length = length(y);

% for saving result
filename = 'tdoa_result_sensor_fixed.mat';

if ~isempty(sensor_position_filename)
    sensor_position = load(sensor_position_filename);
else
    sensor_position = get_sensor_position(sensor_length, uca_radius_meter, sensor_exist_at_uca_center);
end

error_torrieri = nan(x_length, y_length);
for i = 1 : x_length
    for j = 1 : y_length
        
        position_error_torrieri = zeros(1, trial_length);
        target_position = [x(i), y(j)]
        
        if check_target_overlap_sensor(sensor_position, target_position, overlap_criterion_meter)
            continue; % target overlap sensor, so skip
        end
        
        for k = 1 : trial_length
            [position_error_torrieri(k)] = ...
                simulate_tdoa_sensor_fixed(sensor_length, sensor_exist_at_uca_center, ...
                uca_radius_meter, snr_db, target_position, radius_ratio, ...
                sensor_position_filename);
        end
        
        position_error_torrieri;
        mean_torrieri = mean(position_error_torrieri);
        std_torrieri = std(position_error_torrieri);
        
        error_torrieri(i, j) = mean_torrieri;
    end
end
error_torrieri

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

%%
function [] = ...
    plot_spatial_spectrum_2D(spatial_spectrum, delta_search_azimuth, delta_search_elevation, ...
    title_text, figure_position)

if isempty(figure_position)
    figure;
else
    figure('position', figure_position);
end

A = 0:delta_search_azimuth:360 - delta_search_azimuth;
E = 0:delta_search_elevation:90 - delta_search_elevation;
% face lighting with phong make spectrogram more nice
surf(A, E, spatial_spectrum', 'EdgeColor', 'none', 'FaceLighting', 'phong');
% surf(F, T, P', 'EdgeColor', 'none', 'FaceLighting', 'phong');
axis xy; axis tight; colormap(jet); view(0, 90);
ylabel('elevation in deg');
xlabel('azimuth in deg');
if ~isempty(title_text)
    title(title_text);
end

% plot colorbar
h = colorbar;
set(get(h, 'YLabel'), 'String', 'dB');

end






