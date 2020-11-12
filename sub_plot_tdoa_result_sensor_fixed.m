function [] = sub_plot_tdoa_result_sensor_fixed(filename)
% ### called in plot_tdoa_result_sensor_fixed.m

% ### reminder ###
% ### in batch_simulate_tdoa_sensor_fixed.m
% save(filename, 'sensor_position', 'snr_db', 'trial_length', 'x', 'y', 'target_radius_meter', ...
%     'delta_x', 'delta_y', 'error_torrieri', 'cep_mean', 'gdop_mean', 'target_signal_spec', ...
%     'target_position_limit');

load(filename);
% error_torrieri dimension = x_length x y_length
size(error_torrieri);
error_torrieri;

% #### surf dont display last row and column
append_last_row_and_column = 1;

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

% #### after 161109, batch_simulate_tdoa_sensor_fixed.m save 'target_signal_spec', 'target_position_limit'
if exist('target_signal_spec', 'var')
    target_signal_spec
    ndlrb = target_signal_spec(1);
    nprsrb = target_signal_spec(2);
    subframe_length = target_signal_spec(3);
    
    [bw_mhz, fs, nfft, sample_length] = get_bw_from_prs_spec_db(ndlrb, nprsrb, subframe_length);
    bw_mhz, fs, nfft, sample_length
end

figure;
% surf(x, y, error_torrieri', 'EdgeColor', 'none');
surf(x, y, error_torrieri', 'EdgeColor', 'none', 'FaceLighting', 'phong');
axis xy; axis tight; axis equal; colormap(jet); view(0, 90);
xlabel('x axis distance in meter');
ylabel('y axis distance in meter');
title_str = ...
    sprintf('[location error] snr = %d db, trial = %d, target step = [%d, %d] meter', ...
    snr_db, trial_length, delta_x, delta_y);
title(title_str, 'fontweight', 'bold');
% hold on;
annotate_sensor_text(x, y, sensor_position);

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
    sprintf('[cep] snr = %d db, trial = %d, target step = [%d, %d] meter', ...
    snr_db, trial_length, delta_x, delta_y);
title(title_str, 'fontweight', 'bold');
% shading interp;
annotate_sensor_text(x, y, sensor_position);

% plot colorbar
h = colorbar;
set(get(h, 'YLabel'), 'String', 'meter');

% #########################################################
% #### gdop plot is same as cep plot
% #### you can skip to plot gdop: comment out(161018)
% #########################################################

% figure;
% % surf(x, y, error_torrieri', 'EdgeColor', 'none');
% surf(x, y, gdop_mean', 'EdgeColor', 'none', 'FaceLighting', 'phong');
% axis xy; axis tight; axis equal; colormap(jet); view(0, 90);
% xlabel('x axis distance in meter');
% ylabel('y axis distance in meter');
% title_str = ...
%     sprintf('[gdop] snr = %d db, trial = %d, target step = [%d, %d] meter', ...
%     snr_db, trial_length, delta_x, delta_y);
% title(title_str, 'fontweight', 'bold');
% % shading interp;
% annotate_sensor_text(x, y, sensor_position);
% 
% % plot colorbar
% h = colorbar;
% set(get(h, 'YLabel'), 'string', '');

end

