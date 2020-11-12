function [] = ...
    simulate_tdoa_sensor_fixed(sensor_length, sensor_exist_at_uca_center, uca_radius_meter, snr_db, ...
    target_position, radius_ratio, sensor_position_filename, trial_length, fading_param_filename, ...
    xlsrange, target_signal_spec, plot_torrieri_estimation)
% simulate target location using least square estimator
%
% [input]
% - sensor_length: number of sensor, 3 ~ 7
% - sensor_exist_at_uca_center: boolean. 0 = sensor dont exist at uca center, 1 = sensor exist at uca center
% - uca_radius_meter: uca radius in meter, ignored when sensor_position_filename is not empty
% - snr_db: snr in db
% - target_position: target position
% - radius_ratio: ratio between target distance from [0,0] and farthest sensor distance from [0,0]
%   when radius_ratio <= 1, target is located inside sensor boundary
%   when radius_ratio > 1, target can be located outside sensor boundary
% - sensor_position_filename: if not empty, sensor position in text file(.txt) is used
%   sensor_length, uca_radius_meter is ignored
% - trial_length: trial length
% - fading_param_filename: if not empty, fading(rician) is applied, excel file have fading parameter
%   when fading is applied, it is reasonable that target is located within sensor boundary
%   ###### sensor number in fading param file MUST BE SAME as in sensor position file #########
% - xlsrange: excel file range to be read
%   ###############################################
%   #### MUST CHECK excel file(blue bold cell)
%   ###############################################
% - target_signal_spec: target signal(lte prs) spec, vector = [ndlrb, nprsrb, subframe_length]
%   ndlrb and nprsrb determine signal bandwidth, subframe_length determine correlation length
%   for details, see get_bw_from_prs_spec_db.m
% - plot_torrieri_estimation: boolean. 0 = no plot estimation result, 1 = plot estimation result
%
% [usage]
% simulate_tdoa_sensor_fixed(4, 1, 7e3, 5, [100,-100], 1.5, 'sensor_position.txt', 100, [], [], [15,2,1], 1);
% simulate_tdoa_sensor_fixed(4, 1, 7e3, 5, [100,-100], 1.5, 'sensor_position.txt', 100, 'rician_param.xlsx', 'b3:d7', [15,2,1], 1);
% simulate_tdoa_sensor_fixed(4, 1, 7e3, 5, [100,-100], 1.5, [], 100, [], [], [15,2,1], 1);

% #########################################################################
% ### CAUTION: 
% ### when writing code, dont forget target position NOT be (0,0) 
% #########################################################################

%% Check input parameter

% when 0, use all method: torrieri(least square estimator), linear equation solver, hyperbolic equation solver
% when 1, use only torrieri method
% #########################################################################################
% #### MUST BE SET TO 1:
% #### sub_simulate_tdoa_sensor_fixed ONLY RETURN estimation result using torrieri method
% #########################################################################################
use_only_torrieri_method = 1;

% #################################################################################
% ### this is good for sensor planning in real environment
% #################################################################################
if ~isempty(sensor_position_filename)
    [sensor_position, uca_radius_meter] = get_sensor_position_from_file(sensor_position_filename);
    sensor_length = size(sensor_position, 1);
else
    sensor_position = get_sensor_position(sensor_length, uca_radius_meter, sensor_exist_at_uca_center);
end
% sensor_position dimension = sensor_length x 2

if sensor_length > 7 || sensor_length < 3
    fprintf(2, '##### sensor number: 3 ~ 7\n');
    return;
end

% set target position limit 
target_radius_meter = uca_radius_meter * radius_ratio;

if ~check_target_position_is_good(target_position, target_radius_meter)
    fprintf(2, '###### target position is out of range\n');
    return;
end

if ~isempty(fading_param_filename)
    [rician_param, rician_param_raw] = get_rician_parameter_from_excel_file(fading_param_filename, xlsrange);
    rician_param
    if size(rician_param, 1) ~= sensor_length
        fprintf(2, '##### row length of rician parameter must be same as sensor length\n');
        return;
    end
else
    rician_param = [];
end

%%

% ###################################################################################
% #### only when trial_length is very small, set plot_isgnal, plot_position to 1
% #### otherwise numerous figure is generated
% ###################################################################################
plot_signal = 0; % control plot of tx, rx, correlation
plot_position = 0; % control plot of target, sensor, tdoa curve

ndlrb = target_signal_spec(1);
nprsrb = target_signal_spec(2);
subframe_length = target_signal_spec(3);

[bw_mhz, fs, nfft, sample_length] = get_bw_from_prs_spec_db(ndlrb, nprsrb, subframe_length);
bw_mhz, fs, nfft, sample_length

[tx_signal, fs, nfft] = generate_target_signal_lte_prs(ndlrb, nprsrb, subframe_length);
fs, nfft
% tx_signal_length = length(tx_signal);

% ##### generate fsk 422mhz
% [tx_signal, fs, nfft, bw_mhz] = generate_target_signal_fsk_422mhz(fs);

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
        sub_simulate_tdoa_sensor_fixed(sensor_position, snr_db, target_position, ...
        plot_signal, plot_position, use_only_torrieri_method, target_radius_meter, rician_param, ...
        tx_signal, fs, nfft, bw_mhz);
end

x_est_torrieri;
lambda1;
lambda2;
theta_rad;

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

if plot_torrieri_estimation
    if isempty(rician_param)
        title_text = sprintf('[location] snr = %d db, trial = %d, location error = %d m, excl = %d', ...
            snr_db, trial_length, round(mean_position_error), trial_length - sum(idx));
    else
        title_text = sprintf('[location, rician] snr = %d db, trial = %d, location error = %d m, excl = %d', ...
            snr_db, trial_length, round(mean_position_error), trial_length - sum(idx));
    end
    legend_str = plot_sensor_position_only(sensor_position, target_radius_meter, title_text, ...
        fading_param_filename);
    hold on;
    plot_target_position_and_ellipse(target_position, x_est_torrieri, cep, gdop, lambda1, lambda2, theta_rad);
    
    histogram_bin_length = 100;
    plot_position_error_histogram(position_error_excl, histogram_bin_length, title_text, ...
        fading_param_filename);
end

end

%%
% function [] = plot_position_error_histogram(position_error_excl, bin_length, title_text, ...
%     fading_param_filename)
% 
% error_max = max(position_error_excl);
% 
% x = linspace(0, error_max, bin_length);
% 
% N = histc(position_error_excl, x);
% 
% if ~isempty(fading_param_filename)
%     figure('name', fading_param_filename);
% else
%     figure;
% end
% bar(x, N);
% grid on;
% xlabel('position error in meter');
% ylabel('count');
% 
% title(title_text);
% 
% end

%%
% function [] = plot_target_position_and_ellipse(target_position, ...
%     x_est_torrieri, cep, gdop, lambda1, lambda2, theta_rad)
% 
% trial_length = size(x_est_torrieri, 2);
% 
% plot(target_position(1), target_position(2), 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'r', ...
%     'LineWidth', 2);
% 
% for n = 1 : trial_length
%     hold on;
%     plot(x_est_torrieri(1, n), x_est_torrieri(2, n), 'b+', 'LineWidth', 2);
% end
% 
% % concentration ellipse corresponding to probability Pe: 
% % major axes length = 2 * sqrt(k * lambda1) 
% % minor axes length = 2 * sqrt(k * lambda2)
% % where k = -2 * log(1 - Pe)
% % see eq (60)
% 
% Pe = .5;
% k = -2 * log(1 - Pe);
% ra = 2 * sqrt(k * lambda1);
% rb = 2 * sqrt(k * lambda2);
% my_ellipse(ra, rb, theta_rad, x_est_torrieri(1, :), x_est_torrieri(2, :), 'r');
% 
% end

%%
% function [] = plot_sensor_position_only(sensor_position, target_radius_meter, title_text, ...
%     fading_param_filename)
% % copy and modify from plot_sensor_position.m
% 
% colours = my_hPositioningColors();
% 
% sensor_length = size(sensor_position, 1);
% 
% if ~isempty(fading_param_filename)
%     H = figure('name', fading_param_filename);
% else
%     H = figure;
% end
% hold on;
% legendstr = cell(1, sensor_length);
% 
% % Plot the position of each sensor
% for i = 1 : sensor_length
%     plot(sensor_position(i, 1), sensor_position(i, 2), strcat(colours{i}, 'o'), ...
%         'MarkerSize', 7, 'LineWidth', 2);
%     legendstr{i} = sprintf('sensor%d',i);
% end
% 
% % ####### set axis equal #######
% axis equal;
% 
% % ##### when set max axis, consider uca_radius_meter and radius_ratio
% % ##### target_radius_meter = uca_radius_meter * radius_ratio
% axis([-1 1 -1 1] * target_radius_meter * 1.1);
% grid on;
% legend(legendstr);
% xlabel('X position (m)');
% ylabel('Y position (m)');
% title(title_text);
% 
% end



