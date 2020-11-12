function sensor = randomize_sensor_distance_from_target(sensor)

sensor_length = length(sensor);
phi = zeros(sensor_length, 1);
r = zeros(sensor_length, 1);
for n = 1 : sensor_length
%     sensor{n}.Position;
    [phi(n), r(n)] = cart2pol(sensor{n}.Position(1), sensor{n}.Position(2));
end

r = circshift(r, randi([-sensor_length, sensor_length], 1));

for n = 1 : sensor_length
    [sensor{n}.Position(1), sensor{n}.Position(2)] = pol2cart(phi(n), r(n));
end

end