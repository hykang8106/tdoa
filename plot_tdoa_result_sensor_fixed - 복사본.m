function [] = plot_tdoa_result_sensor_fixed(filename)
%
% [usage]
% plot_tdoa_result_sensor_fixed('tdoa_result_sensor_fixed.mat')
%

% ################ results #############
% (1) tdoa_result_sensor_fixed_1610112009.mat: uca(4) + (0,0)
%     as expected 
% (2) tdoa_result_sensor_fixed_1610112335.mat: uca(4) + (0,0)
%     as expected
% (3) tdoa_result_sensor_fixed_1610140302.mat: random(4)
%     ##### NOT as expected: location error is not similar to cep
%     ##### try to add sensor to (0,0)
% (4) 

% #### surf dont display last row and column
append_last_row_and_column = 1;

% ### reminder ###
% ### in batch_simulate_tdoa_sensor_fixed.m
% save(filename, 'sensor_position', 'snr_db', 'trial_length', 'x', 'y', 'target_radius_meter', ...
%     'delta_x', 'delta_y', 'error_torrieri', 'cep_mean', 'gdop_mean');

load(filename);
% error_torrieri dimension = x_length x y_length
size(error_torrieri);
error_torrieri;

if append_last_row_and_column
    x = [x, x(end) + delta_x];
    y = [y, y(end) + delta_y];
    
    % copy and append last row
    error_torrieri = [error_torrieri; error_torrieri(end, :)];
    % copy and append last column
    error_torrieri = [error_torrieri, error_torrieri(:, end)];
    
    size(error_torrieri);
    
    cep_mean = [cep_mean; cep_mean(end, :)];
    cep_mean = [cep_mean, cep_mean(:, end)];
    
    gdop_mean = [gdop_mean; gdop_mean(end, :)];
    gdop_mean = [gdop_mean, gdop_mean(:, end)];
end

figure;
% surf(x, y, error_torrieri', 'EdgeColor', 'none');
surf(x, y, error_torrieri', 'EdgeColor', 'none', 'FaceLighting', 'phong');
axis xy; axis tight; axis equal; colormap(jet); view(0, 90);
xlabel('x axis distance in meter');
ylabel('y axis distance in meter');
title_str = ...
    sprintf('[location error] snr = %d db, trial = %d, distance step = [%d, %d] meter', ...
    snr_db, trial_length, delta_x, delta_y);
title(title_str, 'fontweight', 'bold');
% shading interp;
% hold on;
annotate_sensor_text(error_torrieri, x, y, sensor_position);
% uistack(h_text, 'top');

% plot colorbar
h = colorbar;
set(get(h, 'YLabel'), 'String', 'meter');

figure;
% surf(x, y, error_torrieri', 'EdgeColor', 'none');
surf(x, y, cep_mean', 'EdgeColor', 'none', 'FaceLighting', 'phong');
axis xy; axis tight; axis equal; colormap(jet); view(0, 90);
xlabel('x axis distance in meter');
ylabel('y axis distance in meter');
title_str = ...
    sprintf('[cep] snr = %d db, trial = %d, distance step = [%d, %d] meter', ...
    snr_db, trial_length, delta_x, delta_y);
title(title_str, 'fontweight', 'bold');
% shading interp;
annotate_sensor_text(error_torrieri, x, y, sensor_position);

% plot colorbar
h = colorbar;
set(get(h, 'YLabel'), 'String', 'meter');

figure;
% surf(x, y, error_torrieri', 'EdgeColor', 'none');
surf(x, y, gdop_mean', 'EdgeColor', 'none', 'FaceLighting', 'phong');
axis xy; axis tight; axis equal; colormap(jet); view(0, 90);
xlabel('x axis distance in meter');
ylabel('y axis distance in meter');
title_str = ...
    sprintf('[gdop] snr = %d db, trial = %d, distance step = [%d, %d] meter', ...
    snr_db, trial_length, delta_x, delta_y);
title(title_str, 'fontweight', 'bold');
% shading interp;
annotate_sensor_text(error_torrieri, x, y, sensor_position);

% plot colorbar
h = colorbar;
set(get(h, 'YLabel'), 'string', '');

% figure;
% contour(x, y, error_torrieri');

end

%%
function [] = annotate_sensor_text(error_torrieri, x, y, sensor_position)

sensor_position;

[I, J] = find(isnan(error_torrieri));
sensor_length = length(I);

for n = 1 : sensor_length
    x_position = x(I(n));
    y_position = y(J(n));
    [sensor_number] = ...
        get_sensor_number_from_sensor_position_and_nan_index(x_position, y_position, sensor_position);
    h = text(x_position, y_position, sprintf('s%d', sensor_number), ...
        'VerticalAlignment', 'bottom', 'color', 'r', 'fontweight', 'bold');
    
    % ##### bring 'sensor' text to front
    % ##### reference: https://kr.mathworks.com/matlabcentral/newsreader/view_thread/292321
    set(h,'erasemode','xor');
    set(h,'erasemode','background');
    
%     uistack(h_text, 'top');
end

end

%%
function [sensor_number] = get_sensor_number_from_sensor_position_and_nan_index(x_pos, y_pos, sensor_pos)

% sensor_pos dimension = sensor_length x 2
sensor_length = size(sensor_pos, 1);

[Y, sensor_number] = min(sqrt(sum((repmat([x_pos, y_pos], sensor_length, 1) - sensor_pos).^2, 2)));

end

%% ### surf example

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
