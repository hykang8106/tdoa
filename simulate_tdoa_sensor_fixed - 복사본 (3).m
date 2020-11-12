function [position_error_torrieri, x_est_torrieri, cep, gdop] = ...
    simulate_tdoa_sensor_fixed(sensor_length, sensor_exist_at_uca_center, uca_radius_meter, snr_db, ...
    target_position, radius_ratio, sensor_position_filename, trial_length, fading_param_filename, ...
    target_signal_spec)
% #################################################
% ### rician fading is supported! (161027)
% ### target signal spec is supported! (161029)
% #################################################
%
% simulate target location using tdoa method
%
% this differ from simulate_tdoa_target_fixed:
% ######################################################################
% ##### target position is given by input or randomly determined
% ##### sensors location are fixed: uca + (0,0)
% ##### only torrieri algorithm(least square estimator) is implemented
% ##### array variable is used instead of cell variable
% ######################################################################
%
% [input]
% - sensor_length: number of sensor, 3 ~ 7
% - sensor_exist_at_uca_center: boolean
% - uca_radius_meter: uca radius in meter
%   ###### overriden when sensor_position_filename is not empty
% - snr_db: snr in db
% - target_position: target position
% - radius_ratio: 
%   ### when radius_ratio <= 1, target is located inside uca
%   ### when radius_ratio > 1, target can be located outside uca
% - sensor_position_filename: if not empty, sensor position in file is used, ascii file, .txt
%   ####### sensor_length, uca_radius_meter IS OVERRIDEN
%   ### MUST USE [], not implemented yet
% - trial_length: 
% - fading_param_filename: if not empty, fading(rician) is applied, excel file have fading parameter
%   when fading is applied, it is reasonable that target is located within sensor boundary
%   ###### sensor number in fading param file MUST BE SAME as in sensor position file #########
% - target_signal_spec: target signal(lte prs) spec, vector = [ndlrb, nprsrb, subframe_length]
%   ndlrb and nprsrb determine target signal bandwidth, subframe_length determine target signal length
%   #### for details, see get_bw_from_prs_spec_db.m
% [usage]
% simulate_tdoa_sensor_fixed(4, 1, 7e3, 5, [1055, -3000], 1, [], 100, [], [15,2,1]);
% simulate_tdoa_sensor_fixed(4, 1, 7e3, 5, [1055, -3000], 1, 'sensor_position.txt', 100, [], [15,2,1]);
% simulate_tdoa_sensor_fixed(4, 1, 7e3, 5, [1055, -3000], 1, 'sensor_position.txt', 100, 'rician_param.xlsx', [15,2,1]);
%

% #########################################################################
% ### CAUTION: 
% ### when writing code, dont forget target position NOT be (0,0) 
% #########################################################################

%% Check input parameter

% if sensor_length > 7 || sensor_length < 3
%     fprintf('sensor number: 3 ~ 7\n');
%     return;
% end
% rf_sensor = cell(1,sensor_length);

% sensor_exist_at_uca_center = 1;

plot_torrieri_estimation = 1;

% when 0, use all method: torrieri(least square estimator), linear equation solver, hyperbolic equation solver
% when 1, use only torrieri method
use_only_torrieri_method = 1;

% ###############################################
% #### MUST CHECK excel file(blue bold cell)
% ###############################################
xlsrange = 'b3:d7'; 

% radius_ratio = 1.0; % when radius_ratio <= 1.0, target is located inside uca
% radius_ratio = 1.5;
% radius_ratio = 2.0; % when radius_ratio > 1.0, target can be located outside uca

% if target_position(1) > target_radius_meter || ...
%         target_position(1) < -target_radius_meter || ...
%         target_position(2) > target_radius_meter || ...
%         target_position(2) < -target_radius_meter
%     fprintf('###### target position is out of range\n');
%     return;
% end

% sensor_position, dimension = sensor_length x 2
% #################################################################################
% ### write code to get sensor position from file which have sensor position
% ### this is good for sensor planning in real environment
% #################################################################################
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

if ~check_target_position_is_good(target_position, target_radius_meter)
    fprintf('###### target position is out of range\n');
    return;
end

if ~isempty(fading_param_filename)
    rician_param = get_rician_parameter_from_excel_file(fading_param_filename, xlsrange);
    if size(rician_param, 1) ~= sensor_length
        fprintf('##### row length of rician parameter must be same as sensor length\n');
        return;
    end
else
    rician_param = [];
end

% % get rf sensor position, randomly located
% for i=1:sensor_length
%     rf_sensor{i}.Position = hPositioningPosition(i-1, sensor_length);
%     rf_sensor{i}.Position;
% end

% if isempty(target_position)
%     % target_position, dimension = 1 x 2
%     target_position = random_target_position(uca_radius_meter, radius_ratio);
% end

%%

plot_signal = 0; % control plot of tx, rx, correlation
plot_position = 0; % control plot of target, sensor, tdoa curve

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
        target_signal_spec);
    
    %     if isempty(rician_param)
    %         [position_error_torrieri(n), x_est_torrieri(:, n), cep(n), gdop(n), ...
    %             lambda1(n), lambda2(n), theta_rad(n)] = ...
    %             sub_simulate_tdoa_sensor_fixed(sensor_position, snr_db, target_position, ...
    %             plot_signal, plot_position, use_only_torrieri_method, target_radius_meter);
    %     else
    %         % fading_sub_simulate_tdoa_sensor_fixed is needed
    %         % because it cant be used in batch_simulate_tdoa_sensor_fixed
    %         % batch_simulate_tdoa_sensor_fixed show position error versus all target position
    %         % some target position is not good to apply fading
    %         [position_error_torrieri(n), x_est_torrieri(:, n), cep(n), gdop(n), ...
    %             lambda1(n), lambda2(n), theta_rad(n)] = ...
    %             fading_sub_simulate_tdoa_sensor_fixed(sensor_position, snr_db, target_position, ...
    %             plot_signal, plot_position, use_only_torrieri_method, target_radius_meter, rician_param);
    %     end
end

x_est_torrieri;
lambda1;
lambda2;
theta_rad;

mean_position_error = mean(position_error_torrieri)

if plot_torrieri_estimation
    if isempty(rician_param)
        title_text = sprintf('[location] snr = %d db, trial = %d, mean location error = %f m', ...
            snr_db, trial_length, mean_position_error);
    else
        title_text = sprintf('[location in rician channel] snr = %d db, trial = %d, mean location error = %f m', ...
            snr_db, trial_length, mean_position_error);
    end
    plot_sensor_position_local(sensor_position, target_radius_meter, title_text);
    hold on;
    plot_position_and_ellipse(target_position, x_est_torrieri, cep, gdop, lambda1, lambda2, theta_rad);
    
    bin_length = 100;
    plot_position_error_histogram(position_error_torrieri, bin_length);
end

end

%%
function [] = plot_position_error_histogram(position_error_torrieri, bin_length)

error_max = max(position_error_torrieri);

x = linspace(0, error_max, bin_length);
% delta_x = x(2) - x(1);
% x = [0 : x_step : error_max];

N = histc(position_error_torrieri, x);

figure;
bar(x, N);
grid on;
% xlim([x(1) - delta_x, x(end) + delta_x]);
xlabel('position error in meter');
ylabel('count');

% % define x axis of plot
% x = [0:delta_search_azimuth:360 - delta_search_azimuth];
% 
% % compute histogram of estimated aoa
% estimated_aoa = reshape(estimated_aoa, 1, prod(size(estimated_aoa)));
% N = histc(estimated_aoa, x);
% 
% figure('position', [160 226 735 420]);
% 
% % plot aoa histogram
% subplot(2,1,1);
% bar(x, N);
% grid on;
% xlim([x(1) x(end)]);
% xlabel('azimuth in deg');
% ylabel('count');

end

%%
function [] = plot_position_and_ellipse(target_position, ...
    x_est_torrieri, cep, gdop, lambda1, lambda2, theta_rad)

trial_length = size(x_est_torrieri, 2);

plot(target_position(1), target_position(2), 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'r', ...
    'LineWidth', 2);

for n = 1 : trial_length
    hold on;
    plot(x_est_torrieri(1, n), x_est_torrieri(2, n), 'b+', 'LineWidth', 2);
end

% concentration ellipse corresponding to probability Pe: 
% major axes length = 2 * sqrt(k * lambda1) 
% minor axes length = 2 * sqrt(k * lambda2)
% where k = -2 * log(1 - Pe)
% see eq (60)

Pe = .5;
k = -2 * log(1 - Pe);
ra = 2 * sqrt(k * lambda1);
rb = 2 * sqrt(k * lambda2);
my_ellipse(ra, rb, theta_rad, x_est_torrieri(1, :), x_est_torrieri(2, :), 'r');

end

%%
function [] = plot_sensor_position_local(sensor_position, target_radius_meter, title_text)
% copy from plot_sensor_position.m

colours = my_hPositioningColors();
%     styles = hStyles();

sensor_length = size(sensor_position, 1);

H = figure;
hold on;
legendstr = cell(1, sensor_length);

% Plot the position of each sensor
for i = 1 : sensor_length
    plot(sensor_position(i, 1), sensor_position(i, 2), ...
        strcat(colours{i}, 'o'), ...
        'MarkerSize', 7, 'LineWidth', 2);
    legendstr{i} = sprintf('sensor%d',i);
end

% % Plot the position of target
% plot(target_position(1), target_position(2), 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'r', ...
%     'LineWidth', 2);
% %     plot(0, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'c', ...
% %         'LineWidth', 2);

% ####### set axis equal #######
axis equal;

% ##### when set max axis, consider uca_radius_meter and radius_ratio
axis([-1 1 -1 1] * target_radius_meter * 1.1);
grid on;
% axis([-11000 11000 -11000 11000]);
% legend([legendstr 'target']);
legend(legendstr);
xlabel('X position (m)');
ylabel('Y position (m)');
title(title_text);

end



