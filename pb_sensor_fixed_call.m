function [] = pb_sensor_fixed_call(varargin)
% ################# external function result: FAILED!!

% ### tried to make this code function(gui_get_sensor_position_dialog_box.m),
% ### but when go to "catch" section, sensor_position is still empty
% ### this cause error
if isempty(S.sensor_position)
    uiwait(msgbox({'empty sensor position', 'dialog box will be opened'},...
        'warning','modal'));
    sensor_position_filename = uigetfile('*.txt','select sensor position file');
    try
        [S.sensor_position, S.uca_radius_meter] = ...
            get_sensor_position_from_file(sensor_position_filename);
    catch
        disp('Unable to Load.  Check Name and Try Again.');
        return;
    end
end

% disable load sensor position file and fading parameter file
set(S.fm,{'enable'},{'off';'off'});

[trial_length, radius_ratio, target_pos_y, target_pos_x, ...
    snr_db, subframe_length, nprsrb, ndlrb] = gui_sensor_fixed_input_from_uicontrol(S);

%             trial_length = str2num(get(S.ed(1),'String')); % trial length
%             % get 'string', 'value' property of edit uicontrol
%             get(S.ed(1), {'string', 'value'});
%             % get output example = '100' [0], ##### what is [0]? [ans] no meaning in edit uicontrol, index begin from 1
%
%             radius_ratio = str2num(get(S.ed(2),'String')); % radius ratio
%
%             target_pos_y = str2num(get(S.ed(3),'String')); % y of target position
%
%             target_pos_x = str2num(get(S.ed(4),'String')); % x of target position
%
%             snr_db = str2num(get(S.ed(5),'String')); % snr in db
%
%             % ca = cell array
%             ca = get(S.pp(1), {'string', 'value'}); % subframe length
%             subframe_length = str2num(ca{1}{ca{2}});
%
%             ca = get(S.pp(2), {'string', 'value'}); % nprsrb
%             nprsrb = str2num(ca{1}{ca{2}});
%
%             ca = get(S.pp(3), {'string', 'value'}); % ndlrb
%             ndlrb = str2num(ca{1}{ca{2}});

[bw_mhz, fs, nfft, sample_length] = get_bw_from_prs_spec_db(ndlrb, nprsrb, subframe_length);

[tx_signal, fs, nfft] = generate_target_signal_lte_prs(ndlrb, nprsrb, subframe_length);

plot_signal = 0; plot_position = 0; use_only_torrieri_method = 1;

target_position = [target_pos_x, target_pos_y];

target_radius_meter = S.uca_radius_meter * radius_ratio;

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
        plot_signal, plot_position, use_only_torrieri_method, target_radius_meter, S.rician_param, ...
        tx_signal, fs, nfft, bw_mhz);
end

% ###### target location beyond target_position_limit is ignored:
% ###### excluded in histogram
% ###### excluded in computing mean position error

% #############################################################################################
% ###### current target_position_limit may be too small
% ###### target_position_limit = uca_radius_meter * 3 is good?
% ###### target_position_limit also exist in batch_simulate_tdoa_bw_corr_snr.m
% #############################################################################################
target_position_limit = target_radius_meter;
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
plot_sensor_position_only(S.sensor_position, target_radius_meter, title_text, ...
    S.fading_param_filename);
hold on;
plot_target_position_and_ellipse(target_position, x_est_torrieri, ...
    cep, gdop, lambda1, lambda2, theta_rad);

histogram_bin_length = 100;
plot_position_error_histogram(position_error_excl, histogram_bin_length, title_text, ...
    S.fading_param_filename);

% enable load sensor position file and fading parameter file
set(S.fm,{'enable'},{'on';'on'});

end