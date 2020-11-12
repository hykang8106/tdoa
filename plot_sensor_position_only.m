function [legendstr] = plot_sensor_position_only(sensor_position, target_radius_meter, title_text, ...
    fading_param_filename)
% copy and modify from plot_sensor_position.m

colours = my_hPositioningColors();

sensor_length = size(sensor_position, 1);

if ~isempty(fading_param_filename)
    H = figure('name', fading_param_filename);
else
    H = figure;
end
hold on;
legendstr = cell(1, sensor_length);

% Plot the position of each sensor
for i = 1 : sensor_length
    plot(sensor_position(i, 1), sensor_position(i, 2), strcat(colours{i}, 'o'), ...
        'MarkerSize', 7, 'LineWidth', 2);
    legendstr{i} = sprintf('sensor%d',i);
end

% ####### set axis equal #######
axis equal;

% ##### when set max axis, consider uca_radius_meter and radius_ratio
% ##### target_radius_meter = uca_radius_meter * radius_ratio
axis([-1 1 -1 1] * target_radius_meter * 1.1);
grid on;
% legend(legendstr);
% ###################
% Starting in R2017a, the legend automatically updates when you add or remove data series from the axes. 
% If you do not want the legend to automatically update, set the AutoUpdate property of the legend to 'off'.
% ####################
% legend(legendstr, 'AutoUpdate', 'off');
xlabel('X position (m)');
ylabel('Y position (m)');
title(title_text);

end

