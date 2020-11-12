function [sensor_position] = shrink_sensor_position(sensor_position_filename, shrink_ratio)
% shrink or expand sensor position from file, and save sensor position to new file
%
% [input]
% - sensor_position_filename:
% - shrink_ratio: when shrink_ratio < 1, shrink sensor position. when shrink_ratio > 1, expand sensor position
% [usage]
% shrink_sensor_position('sensor_position.txt', .5);

sp = load(sensor_position_filename);

[theta, r] = cart2pol(sp(:, 1), sp(:, 2));

r = r * shrink_ratio;

[X, Y] = pol2cart(theta, r);

sensor_position = [X, Y];

[p, name, e] = fileparts(sensor_position_filename);
new_filename = [name, '_', num2str(shrink_ratio), e]

save(new_filename, 'sensor_position', '-ascii');

end