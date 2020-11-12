function plot_xcorrelation(corr, fs, ref_sensor_number)

colors = my_hPositioningColors();
figure;
hold on;

[sample_length, sensor_length] = size(corr);

% Plot correlation for each eNodeB
t = (0 : sample_length - 1) / fs;
legendstr = cell(1, sensor_length);
for i = 1 : sensor_length
    plot(t, abs(corr(:, i)), colors{i}, 'LineWidth', 2);
    legendstr{i} = sprintf('sensor%d', i);
end
grid on;
xlim([t(1) t(end)]);
legend(legendstr);
xlabel('Time (s)');
ylabel('Absolute value');
title(sprintf('Receiver correlations, reference sensor = %d', ref_sensor_number));

% Plot markers are correlation peaks
for i = 1 : sensor_length
    c = abs(corr(:, i));
    j = find(c == max(c), 1);
    plot(t(j), c(j), strcat(colors{i}, 'o'), 'MarkerSize', 5, ...
        'MarkerFaceColor', colors{i}, 'LineWidth', 2);
    legendstr{i} = sprintf('sensor%d', i);
end

end