function plot_rx_signal(rx, fs, snr_db)

colors = my_hPositioningColors();
figure
hold on;

[sample_length, sensor_length] = size(rx);
t = (0 : sample_length - 1) / fs;  % Calculate time of samples

% For each eNodeB plot the received waveform
legendstr = cell(1, sensor_length);
for i = 1 : sensor_length
    plot(t, abs(rx(:, i)), colors{i}, 'LineWidth', 2);
    legendstr{i} = sprintf('sensor%d',i);
end

% legendstr = cell(1, numel(rx));
% for i = 1:numel(rx)
%     plot(t, abs(rx{i}), colors{i}, 'LineWidth', 2);
%     legendstr{i} = sprintf('sensor%d',i);
% end

grid on;
xlim([t(1) t(end)]);
legend(legendstr);
xlabel('Time (s)');
ylabel('Absolute value');
title(sprintf('Received waveforms at each sensor location, snr = %d db', snr_db));

end
