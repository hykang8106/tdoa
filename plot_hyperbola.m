function [] = plot_hyperbola(sensor1_position, delta_distance, axis_equal)
% plot hyperbola
%
% [input]
% - sensor1_position: located on negative x-axis, [-7, 0]
% - delta_distance: distance difference between target from sensor 1 and sensor 2
%   positive = target is near sensor 2
%   negative = target is near sensor 1
% - axis_equal: set x, y axis ratio to equal, 0 or 1
% [usage]
% plot_hyperbola([-7,0], 2, 0)

% see https://en.wikipedia.org/wiki/Hyperbola

% also located on positive x-axis, opposite to sensor 1
sensor2_position = sensor1_position * -1;

% hyperbola parameter
c = sensor1_position(1); % hyperbola focus, x component of sensor 2 position
a = delta_distance/2; % hyperbola vertex, delta distance is equal to distance between two hyperbola vertex
b = sqrt(c^2 - a^2);

% hyperbola equation: (x/a)^2 - (y/b)^2 = 1
% hyperbola parametric equation
t = (-3:.01:3)';    % parameter
x = a * cosh(t);
y = b * sinh(t);

% figure;
% plot(t, [x,y]);

x_max = max(abs(x));
% x_min = min(x);
y_max = max(abs(y));
% y_min = min(y);

figure;
plot(x, y);
axis([-x_max x_max -y_max y_max]);
if axis_equal
    axis equal;
end
grid on;

hold on;
plot(sensor1_position(1), sensor1_position(2), 'sb', ...
    'MarkerSize', 7, 'LineWidth', 2);
plot(sensor2_position(1), sensor2_position(2), 'sm', ...
    'MarkerSize', 7, 'LineWidth', 2);
legend('tdoa line', 'sensor 1', 'sensor 2');

end

