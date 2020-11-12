function [] = plot_tdoa_result(filename)
%
% [usage]
% plot_tdoa_result('tdoa_result.mat')

% ### reminder of file contents
% save(filename, 'sensor_length', 'snr_db', 'trial_length', ...
%     'torrieri_mean', 'torrieri_std', ...
%     'hyperbolic_mean', 'hyperbolic_std', ...
%     'linear_mean', 'linear_std');

load(filename);

s_len = length(sensor_length);
legend_str = cell(s_len, 1);
for n = 1 : s_len
    legend_str{n} = sprintf('%d sensors', sensor_length(n));
end
figure('Position', [188 430 974 420]);
subplot(1,2,1); plot(snr_db, torrieri_mean, 's-'); grid on;
xlabel('snr in db'); 
xlim([snr_db(1), snr_db(end)]); 
ylabel('meter');
legend(legend_str);
title(sprintf('mean of position error: method = torrieri, trial = %d', trial_length));
subplot(1,2,2); plot(snr_db, torrieri_std, 's-'); grid on;
xlabel('snr in db'); 
xlim([snr_db(1), snr_db(end)]); 
ylabel('meter');
legend(legend_str);
title(sprintf('std of position error: method = torrieri, trial = %d', trial_length));

figure('Position', [188 430 974 420]);
subplot(1,2,1); plot(snr_db, hyperbolic_mean, 's-'); grid on;
xlabel('snr in db'); xlim([snr_db(1), snr_db(end)]); ylabel('meter');
legend('3 sensors');
title(sprintf('mean of position error: method = hyperbolic, trial = %d', trial_length));
subplot(1,2,2); plot(snr_db, hyperbolic_std, 's-'); grid on;
xlabel('snr in db'); xlim([snr_db(1), snr_db(end)]); ylabel('meter');
legend('3 sensors');
title(sprintf('std of position error: method = hyperbolic, trial = %d', trial_length));

I = find(sensor_length >= 4);
s_len = length(I);
legend_str = cell(s_len, 1);
for n = 1 : s_len
    legend_str{n} = sprintf('%d sensors', sensor_length(I(n)));
end
figure('Position', [188 430 974 420]);
subplot(1,2,1); plot(snr_db, linear_mean, 's-'); grid on;
xlabel('snr in db'); xlim([snr_db(1), snr_db(end)]); ylabel('meter');
legend(legend_str);
title(sprintf('mean of position error: method = linear, trial = %d', trial_length));
subplot(1,2,2); plot(snr_db, linear_std, 's-'); grid on;
xlabel('snr in db'); xlim([snr_db(1), snr_db(end)]); ylabel('meter');
legend(legend_str);
title(sprintf('std of position error: method = linear, trial = %d', trial_length));

figure;
idx_sensor_length = 1;
errorbar(snr_db, torrieri_mean(:, idx_sensor_length), torrieri_std(:, idx_sensor_length)/2, 's-'); 
grid on;
xlabel('snr in db'); ylabel('meter');
title(sprintf('mean and std of position error: method = torrieri, trial = %d', trial_length));
legend(sprintf('%d sensors', sensor_length(idx_sensor_length)));

end
