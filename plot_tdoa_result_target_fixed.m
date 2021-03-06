function [] = plot_tdoa_result_target_fixed(filename, sensor_number_torrieri)
% plot result of batch_simulate_tdoa_target_fixed.m
%
% [input]
% - filename: if empty, dialog box is displayed for file selection
%   tdoa result file name which include 'tdoa_result_target_fixed' string
% - sensor_number_torrieri: specific sensor number for least square estimator
%   mean and std is plotted together in figure using errorbar
%
% [usage]
% plot_tdoa_result_target_fixed('', 5)

% ### reminder of file contents
% save(filename, 'sensor_length', 'snr_db', 'trial_length', ...
%     'torrieri_mean', 'torrieri_std', ...
%     'hyperbolic_mean', 'hyperbolic_std', ...
%     'linear_mean', 'linear_std');

if isempty(filename)
    filterspec = 'tdoa_result_target_fixed_*.mat';
    [filename, pathname, filterindex] = uigetfile(filterspec);
    if ~filename
        fprintf(2, '######## file selection canceled\n');
        return;
    else
        fprintf('filename = %s\n', filename);
    end
end

sub_plot_tdoa_result_target_fixed(filename, sensor_number_torrieri);

% load(filename);
% 
% filename, sensor_length, snr_db, trial_length
% 
% s_len = length(sensor_length);
% legend_str = cell(s_len, 1);
% for n = 1 : s_len
%     legend_str{n} = sprintf('%d sensors', sensor_length(n));
% end
% 
% if ~isempty(torrieri_mean)
%     figure('Position', [188 430 974 420]);
%     subplot(1,2,1); plot(snr_db, torrieri_mean, 's-'); grid on;
%     xlabel('snr in db');
%     xlim([snr_db(1), snr_db(end)]);
%     ylabel('meter');
%     legend(legend_str);
%     title(sprintf('mean of location error: method = torrieri, trial = %d', trial_length));
%     subplot(1,2,2); plot(snr_db, torrieri_std, 's-'); grid on;
%     xlabel('snr in db');
%     xlim([snr_db(1), snr_db(end)]);
%     ylabel('meter');
%     legend(legend_str);
%     title(sprintf('std of location error: method = torrieri, trial = %d', trial_length));
% end
% 
% if ~isempty(hyperbolic_mean)
%     figure('Position', [188 430 974 420]);
%     subplot(1,2,1); plot(snr_db, hyperbolic_mean, 's-'); grid on;
%     xlabel('snr in db'); xlim([snr_db(1), snr_db(end)]); ylabel('meter');
%     legend('3 sensors');
%     title(sprintf('mean of location error: method = hyperbolic, trial = %d', trial_length));
%     subplot(1,2,2); plot(snr_db, hyperbolic_std, 's-'); grid on;
%     xlabel('snr in db'); xlim([snr_db(1), snr_db(end)]); ylabel('meter');
%     legend('3 sensors');
%     title(sprintf('std of location error: method = hyperbolic, trial = %d', trial_length));
% end
% 
% I = find(sensor_length >= 4);
% s_len = length(I);
% legend_str = cell(s_len, 1);
% for n = 1 : s_len
%     legend_str{n} = sprintf('%d sensors', sensor_length(I(n)));
% end
% 
% if ~isempty(linear_mean)
%     figure('Position', [188 430 974 420]);
%     subplot(1,2,1); plot(snr_db, linear_mean, 's-'); grid on;
%     xlabel('snr in db'); xlim([snr_db(1), snr_db(end)]); ylabel('meter');
%     legend(legend_str);
%     title(sprintf('mean of location error: method = linear, trial = %d', trial_length));
%     subplot(1,2,2); plot(snr_db, linear_std, 's-'); grid on;
%     xlabel('snr in db'); xlim([snr_db(1), snr_db(end)]); ylabel('meter');
%     legend(legend_str);
%     title(sprintf('std of location error: method = linear, trial = %d', trial_length));
% end
% 
% I = find(sensor_length == sensor_number_torrieri);
% if isempty(I)
%     fprintf(2, '#### sensor number for torrieri NOT exist\n');
% end
% 
% if ~isempty(torrieri_mean) && ~isempty(I)
%     figure;
%     errorbar(snr_db, torrieri_mean(:, I), torrieri_std(:, I)/2, 's-');
%     grid on;
%     xlabel('snr in db'); ylabel('meter');
%     title(sprintf('mean and std of location error: method = torrieri, trial = %d, senor = %d', ...
%         trial_length, sensor_number_torrieri));
%     legend(sprintf('%d sensors', sensor_length(I)));
% end

end
