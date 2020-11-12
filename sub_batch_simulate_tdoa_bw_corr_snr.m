function [filename] = sub_batch_simulate_tdoa_bw_corr_snr(sensor_position, uca_radius_meter, ...
    target_position, radius_ratio, trial_length, rician_param, ndlrb, snr_vec, nprsrb_vec, nsubframe_vec)
% subroutine

% when 0, use all method: torrieri(least square estimator), linear equation solver, hyperbolic equation solver
% when 1, use only torrieri method
use_only_torrieri_method = 1;

target_radius_meter = uca_radius_meter * radius_ratio;

% ###### target location beyond target_position_limit is ignored: 
% ###### excluded in histogram
% ###### excluded in computing mean position error

% #############################################################################################
% ###### current target_position_limit is too small, 
% ###### target_position_limit = uca_radius_meter * 3 is good?
% ###### target_position_limit also exist in simulate_tdoa_sensor_fixed.m
% #############################################################################################
target_position_limit = target_radius_meter;

if ~check_target_position_is_good(target_position, target_radius_meter)
    fprintf(2, '###### target position is out of range\n');
    return;
end

% ##############################################################################
% #### dont set to 1 to prevent numerous figure from generating
% ##############################################################################
plot_signal = 0; % control plot of tx, rx, correlation
plot_position = 0; % control plot of target, sensor, tdoa curve

snr_len = length(snr_vec);
nprsrb_len = length(nprsrb_vec);
nsubframe_len = length(nsubframe_vec);

% ##### tip: for dimension rearrange, use "permute" function
position_error_array = zeros(nprsrb_len, snr_len, nsubframe_len);

error_cell_array = cell(nprsrb_len, snr_len, nsubframe_len);

h = waitbar(0, 'Please wait...', 'name', 'batch program progress', 'windowstyle', 'modal');
steps = nprsrb_len * snr_len * nsubframe_len;

u = tic;

for i = 1 : nsubframe_len
    for j = 1 : snr_len
        for k = 1 : nprsrb_len
            
            snr_db = snr_vec(j);
            
            nprsrb = nprsrb_vec(k);
            subframe_length = nsubframe_vec(i);
            
            [bw_mhz, fs, nfft, sample_length] = get_bw_from_prs_spec_db(ndlrb, nprsrb, subframe_length);
            bw_mhz, fs, nfft, sample_length
            
            [tx_signal, fs, nfft] = generate_target_signal_lte_prs(ndlrb, nprsrb, subframe_length);
            fs, nfft
            
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
            
            idx = position_error_torrieri <= target_position_limit;
            position_error_excl = position_error_torrieri(idx);
            
            mean_position_error = mean(position_error_excl);
            
            position_error_array(k, j, i) = mean_position_error;
            
            error_cell_array{k, j, i} = position_error_excl;
            
            step = snr_len * nprsrb_len * (i - 1) + nprsrb_len * (j - 1) + k;
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

% for saving result
filename = 'tdoa_result_bw_corr_snr.mat';

% #### append date string to filename
[p, name, e] = fileparts(filename);
filename = [name, '_', datestr(now, 'yymmddHHMM'), e];
save(filename, 'sensor_position', 'target_position', 'trial_length', 'rician_param', ...
    'ndlrb', 'snr_vec', 'nprsrb_vec', 'nsubframe_vec', 'position_error_array', 'target_position_limit', ...
    'error_cell_array', 'radius_ratio');

end

