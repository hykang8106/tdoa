function [] = ...
    batch_simulate_tdoa_sensor_fixed(sensor_length, sensor_exist_at_uca_center, uca_radius_meter, ...
    snr_db, radius_ratio, sensor_position_filename, trial_length, delta_distance, target_signal_spec)
% batch version of simulate_tdoa_sensor_fixed.m
% #########################################################################################
% ### how location error is varied according to target position and sensor length/shape?
% #########################################################################################
%
% [input]
% - sensor_length: number of sensor, 3 ~ 7
% - sensor_exist_at_uca_center: boolean. 0 = sensor dont exist at uca center, 1 = sensor exist at uca center
% - uca_radius_meter: uca radius in meter, ignored when sensor_position_filename is not empty
% - snr_db: snr in db
% - radius_ratio: ratio between target distance from [0,0] and farthest sensor distance from [0,0]
%   when radius_ratio <= 1, target is located inside sensor boundary
%   when radius_ratio > 1, target can be located outside sensor boundary
% - sensor_position_filename: if not empty, sensor position in text file(.txt) is used
%   sensor_length, uca_radius_meter, sensor_exist_at_uca_center is ignored
% - trial_length: trial length
% - delta_distance: target position delta distance. vector = [delta_x, delta_y]
% - target_signal_spec: target signal(lte prs) spec, vector = [ndlrb, nprsrb, subframe_length]
%   ndlrb and nprsrb determine signal bandwidth, subframe_length determine signal sample length
%   for details, see get_bw_from_prs_spec_db.m
%
% [usage]
% batch_simulate_tdoa_sensor_fixed(5, 1, 5e3, 5, 1.5, 'sensor_position.txt', 300, [100,100], [15,2,1])
% batch_simulate_tdoa_sensor_fixed(5, 1, 5e3, 5, 1.5, [], 100, [100,100], [15,2,1])

% ############################
% #### MUST be empty!!!
% ############################
% fading cant be applied in batch_simulate_tdoa_sensor_fixed
% because batch_simulate_tdoa_sensor_fixed show position error versus all target position
% target position near sensor is not good to apply fading
rician_param = [];

% ##############################################################################
% #### dont set to 1 to prevent numerous figure from generating
% ##############################################################################
plot_signal = 0;
plot_position = 0;

use_only_torrieri_method = 1;

% used to determine whether target is overlapped with sensor
overlap_criterion_meter = 1;
% overlap_criterion_meter = 10;
% overlap_criterion_meter = physconst('LightSpeed') / fs / 2;

delta_x = delta_distance(1);
delta_y = delta_distance(2);

% for saving result
filename = 'tdoa_result_sensor_fixed.mat';

if ~isempty(sensor_position_filename)
    [sensor_position, uca_radius_meter] = get_sensor_position_from_file(sensor_position_filename);
    sensor_length = size(sensor_position, 1);
else
    sensor_position = get_sensor_position(sensor_length, uca_radius_meter, sensor_exist_at_uca_center);
end

if sensor_length > 7 || sensor_length < 3
    fprintf(2, '#### sensor number: 3 ~ 7\n');
    return;
end

target_radius_meter = uca_radius_meter * radius_ratio;
% ##### target_position_limit: differ from batch_simulate_tdoa_bw_corr_snr.m
target_position_limit = uca_radius_meter * 3;

x_range = [-1, 1] * target_radius_meter;
y_range = [-1, 1] * target_radius_meter;

x = x_range(1) : delta_x : x_range(end);
x_length = length(x);
y = y_range(1) : delta_y : y_range(end);
y_length = length(y);

ndlrb = target_signal_spec(1);
nprsrb = target_signal_spec(2);
subframe_length = target_signal_spec(3);

[bw_mhz, fs, nfft, sample_length] = get_bw_from_prs_spec_db(ndlrb, nprsrb, subframe_length);

[tx_signal, fs, nfft] = generate_target_signal_lte_prs(ndlrb, nprsrb, subframe_length);

error_torrieri = nan(x_length, y_length);
cep_mean = nan(x_length, y_length);
gdop_mean = nan(x_length, y_length);

h = waitbar(0, 'Please wait...', 'name', 'batch program progress', 'windowstyle', 'modal');
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
            fprintf(2, '###### target position is out of range\n');
            return;
        end
        
        if check_target_overlap_sensor(sensor_position, target_position, overlap_criterion_meter)
            % target overlap sensor, so skip
            continue; 
        end
        
        for k = 1 : trial_length
            [position_error_torrieri(k), x_est_torrieri(:, k), cep(k), gdop(k), ...
                lambda1, lambda2, theta_rad] = ...
                sub_simulate_tdoa_sensor_fixed(sensor_position, snr_db, target_position, ...
                plot_signal, plot_position, use_only_torrieri_method, target_radius_meter, rician_param, ...
                tx_signal, fs, nfft, bw_mhz);
        end
        
        idx = position_error_torrieri <= target_position_limit;
        position_error_excl = position_error_torrieri(idx);
        
        mean_torrieri = mean(position_error_excl);
        std_torrieri = std(position_error_excl);
        
        error_torrieri(i, j) = mean_torrieri;
        cep_mean(i, j) = mean(cep(idx));
        gdop_mean(i, j) = mean(gdop(idx));
        
        step = y_length * (i - 1) + j;
        w = toc(u);
        % reference: 
        % http://stackoverflow.com/questions/12210583/is-there-a-matlab-function-to-convert-elapsed-seconds-to-hhmmss-format
        z = fix(mod(w, [0, 3600, 60])./[3600, 60, 1]);
        waitbar(step / steps, h, sprintf('%.2f %%, elapsed time = %d : %02d : %02d', ...
            step / steps * 100, z(1), z(2), z(3)));
        
    end
end

close(h);

error_torrieri;
cep_mean;
gdop_mean;

% #### append date string to filename
[p, name, e] = fileparts(filename);
filename = [name, '_', datestr(now, 'yymmddHHMM'), e];
save(filename, 'sensor_position', 'snr_db', 'trial_length', 'x', 'y', 'target_radius_meter', ...
    'delta_x', 'delta_y', 'error_torrieri', 'cep_mean', 'gdop_mean', 'target_signal_spec', ...
    'target_position_limit');

end

