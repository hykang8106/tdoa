function [sensor_position, uca_radius_meter] = get_sensor_position_from_file(sensor_position_filename)

% sensor_position dimension = sensor_length x 2
sensor_position = load(sensor_position_filename);

uca_radius_meter = max(sqrt(sum(sensor_position.^2, 2)));

% sensor_length = size(sensor_position, 1);
% 
% R = zeros(sensor_length, 1);
% 
% for n = 1 : sensor_length
%     R(n) = sqrt(sum(sensor_position(n, :).^2));
% end

end
