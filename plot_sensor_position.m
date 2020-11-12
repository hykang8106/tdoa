function [H] = plot_sensor_position(sensor_position, target_position, target_radius_meter)

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

% Plot the position of target
plot(target_position(1), target_position(2), 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'r', ...
    'LineWidth', 2);
%     plot(0, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'c', ...
%         'LineWidth', 2);

% ##### when set max axis, consider uca_radius_meter and radius_ratio
axis([-1 1 -1 1] * target_radius_meter * 1.1);
% axis([-11000 11000 -11000 11000]);
legend([legendstr 'target']);
xlabel('X position (m)');
ylabel('Y position (m)');
title('Positions');

end