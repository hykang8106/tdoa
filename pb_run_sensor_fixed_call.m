function [] = pb_run_sensor_fixed_call(varargin)

hObject = varargin{3};
S = guidata(hObject);

% ### tried to make this code function(gui_get_sensor_position_dialog_box.m),
% ### but when go to "catch" section, sensor_position is still empty
% ### this cause error
if isempty(S.sensor_position)
    uiwait(msgbox({'empty sensor position: load sensor postion'},...
        'warning','modal'));
    return;
end

if isempty(S.tx_signal)
    uiwait(msgbox({'empty tx signal: generate tx signal'},...
        'warning','modal'));
    return;
end

% disable all uimenu
set(S.fm,{'enable'},repmat({'off'},length(S.fm),1));

%             set(S.fh(2), 'windowstyle', 'modal');

[trial_length, radius_ratio, target_pos_y, target_pos_x, snr_db] = ...
    gui_sensor_fixed_input_from_uicontrol(S);

plot_signal = 0; plot_position = 0; use_only_torrieri_method = 1;

target_position = [target_pos_x, target_pos_y];

S.target_radius_meter = S.uca_radius_meter * radius_ratio;

position_error_torrieri = zeros(1, trial_length);
cep = zeros(1, trial_length);
gdop = zeros(1, trial_length);
lambda1 = zeros(1, trial_length);
lambda2 = zeros(1, trial_length);
theta_rad = zeros(1, trial_length);
% ###### for computing position estimation covariance
x_est_torrieri = zeros(2, trial_length);

for n = 1 : trial_length
    [position_error_torrieri(n), x_est_torrieri(:, n), cep(n), gdop(n), ...
        lambda1(n), lambda2(n), theta_rad(n)] = ...
        sub_simulate_tdoa_sensor_fixed(S.sensor_position, snr_db, target_position, ...
        plot_signal, plot_position, use_only_torrieri_method, S.target_radius_meter, S.rician_param, ...
        S.tx_signal, S.fs, S.nfft, S.bw_mhz);
end

% ###### target location beyond target_position_limit is ignored:
% ###### excluded in histogram
% ###### excluded in computing mean position error

% #############################################################################################
% ###### current target_position_limit may be too small
% ###### target_position_limit = uca_radius_meter * 3 is good?
% ###### target_position_limit also exist in batch_simulate_tdoa_bw_corr_snr.m
% #############################################################################################
target_position_limit = S.target_radius_meter;
idx = position_error_torrieri <= target_position_limit;
position_error_excl = position_error_torrieri(idx);

mean_position_error = mean(position_error_excl);

if isempty(S.rician_param)
    title_text = sprintf('[location] snr = %d db, trial = %d, location error = %d m, excl = %d', ...
        snr_db, trial_length, round(mean_position_error), trial_length - sum(idx));
else
    title_text = sprintf('[location, rician] snr = %d db, trial = %d, location error = %d m, excl = %d', ...
        snr_db, trial_length, round(mean_position_error), trial_length - sum(idx));
end
plot_sensor_position_only(S.sensor_position, S.target_radius_meter, title_text, ...
    S.fading_param_filename);
hold on;
plot_target_position_and_ellipse(target_position, x_est_torrieri, ...
    cep, gdop, lambda1, lambda2, theta_rad);

histogram_bin_length = 100;
plot_position_error_histogram(position_error_excl, histogram_bin_length, title_text, ...
    S.fading_param_filename);

%             set(S.fh(2), 'windowstyle', 'normal');

% enable all uimenu
set(S.fm,{'enable'},repmat({'on'},length(S.fm),1));

guidata(hObject, S);

end

