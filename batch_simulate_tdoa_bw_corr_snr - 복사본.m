function [] = batch_simulate_tdoa_bw_corr_snr(sensor_length, sensor_exist_at_uca_center, uca_radius_meter, ...
    target_position, radius_ratio, sensor_position_filename, trial_length, fading_param_filename, ...
    xlsrange, ndlrb, snr_vec, nprsrb_vec, nsubframe_vec)
% batch version of simulate_tdoa_sensor_fixed.m
% i want to see how position error is varied according to target signal bw, correlation length, and snr
% target and sensor is fixed
%
% [input]
% - sensor_length: number of sensor, 3 ~ 7
% - sensor_exist_at_uca_center: boolean
% - uca_radius_meter: uca radius in meter
%   ###### overriden when sensor_position_filename is not empty
% - target_position: target position. vector, i.e. [0,0]
% - radius_ratio: 
%   ### when radius_ratio <= 1, target is located inside uca
%   ### when radius_ratio > 1, target can be located outside uca
% - sensor_position_filename: if not empty, sensor position in file is used, ascii file, .txt
%   ####### sensor_length, uca_radius_meter, sensor_exist_at_uca_center IS OVERRIDDEN
% - trial_length: 
% - fading_param_filename: if not empty, fading(rician) is applied, excel file have fading parameter
%   when fading is applied, it is reasonable that target is located within sensor boundary
%   ###### sensor number in fading param file MUST BE SAME as in sensor position file #########
% - xlsrange: range to be read in excel file
%   ###############################################
%   #### MUST CHECK excel file(blue bold cell)
%   ###############################################
% - ndlrb: number of downlink resource block in lte prs
% - snr_vec: snr vector in db
% - nprsrb_vec: nprsrb(number of prs resource block) vector, determine target signal bw
%   #### for details, see get_bw_from_prs_spec_db.m
% - nsubframe_vec: nsubframe(number of subframe) vector, determine correlation length
%   #### for details, see get_bw_from_prs_spec_db.m
%
% [usage]
% batch_simulate_tdoa_bw_corr_snr(5, 0, 4e3, [0,0], 1.5, 'sensor_position.txt', 100, [], [], 15, [0,5,10], [1:2:15], [1,2])
% batch_simulate_tdoa_bw_corr_snr(5, 0, 4e3, [0,0], 1.5, 'sensor_position.txt', 100, 'rician_param.xlsx', 'b3:d7', 15, [0,5,10], [1:2:15], [1,2])
%
% ###################################################################################################
% ### test results of no fading: USEFUL!
% ###
% ### test results of rician fading: NOT USEFUL! revise another scheme
% ### to save run time, put routine to read rician parameter from excel file to outside of for loop
% ###################################################################################################

% for saving result
filename = 'tdoa_result_bw_corr_snr.mat';

snr_len = length(snr_vec);
nprsrb_len = length(nprsrb_vec);
nsubframe_len = length(nsubframe_vec);

% ######### tip: for dimension rearrange, use "permute" function
position_error_array = zeros(nprsrb_len, snr_len, nsubframe_len);

% ############################################################
% ### if set to 1, MUST MAKE LOOP VERY SMALL
% ### otherwise large number of matlab figure is generated
% ############################################################
plot_torrieri_estimation = 0;

h = waitbar(0, 'Please wait...', 'name', 'batch program progress');
% hw = findobj(h,'Type','Patch');
% set(hw,'EdgeColor',[0 1 0],'FaceColor',[0 1 0]) % changes the color to green
steps = nprsrb_len * snr_len * nsubframe_len;

u = tic;

for i = 1 : nsubframe_len
    for j = 1 : snr_len
        for k = 1 : nprsrb_len
            
            snr_db = snr_vec(j);
            % target_signal_spec: target signal(lte prs) spec, vector = [ndlrb, nprsrb, subframe_length]
            target_signal_spec = [ndlrb, nprsrb_vec(k), nsubframe_vec(i)];
            
            [mean_position_error, position_error_torrieri, x_est_torrieri, cep, gdop, ...
                sensor_position, rician_param] = ...
                simulate_tdoa_sensor_fixed(sensor_length, sensor_exist_at_uca_center, uca_radius_meter, ...
                snr_db, ...
                target_position, radius_ratio, sensor_position_filename, trial_length, ...
                fading_param_filename, xlsrange, ...
                target_signal_spec, plot_torrieri_estimation);

            % #####################################################################################################
            % ########## consider to use sub_simulate_tdoa_sensor_fixed function to save program run time
            % #####################################################################################################
%             [position_error_torrieri(n), x_est_torrieri(:, n), cep(n), gdop(n), ...
%                 lambda1(n), lambda2(n), theta_rad(n)] = ...
%                 sub_simulate_tdoa_sensor_fixed(sensor_position, snr_db, target_position, ...
%                 plot_signal, plot_position, use_only_torrieri_method, target_radius_meter, rician_param, ...
%                 tx_signal, fs, nfft, bw_mhz);
            
            position_error_array(k, j, i) = mean_position_error;
            
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

% #### append date string to filename
[p, name, e] = fileparts(filename);
filename = [name, '_', datestr(now, 'yymmddHHMM'), e];
save(filename, 'sensor_position', 'target_position', 'trial_length', 'rician_param', ...
    'ndlrb', 'snr_vec', 'nprsrb_vec', 'nsubframe_vec', 'position_error_array');
% save(filename, 'sensor_position', 'target_position', 'trial_length', 'fading_param_filename', ...
%     'ndlrb', 'snr_vec', 'nprsrb_vec', 'nsubframe_vec', 'position_error_array');

end


