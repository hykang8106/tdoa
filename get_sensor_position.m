function [sensor_position] = ...
    get_sensor_position(sensor_length, uca_radius_meter, sensor_exist_at_uca_center)

if sensor_exist_at_uca_center
    uca_length = sensor_length - 1;
else
    uca_length = sensor_length;
end

theta = linspace(0, 2 * pi, uca_length + 1)';

x = uca_radius_meter * cos(theta(1 : end - 1));
y = uca_radius_meter * sin(theta(1 : end - 1));

if sensor_exist_at_uca_center
    % add sensor at (0, 0)
    sensor_position = [0, 0; x, y];
else
    sensor_position = [x, y];
end
sensor_position;

% % add sensor at (0, 0)
% sensor_position = [0, 0; x, y];
% sensor_position

end