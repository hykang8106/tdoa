function [] = pah(sensor1_position, sensor2_position, delta_distance)
% plot arbitrary hyperbola
% generalized version of plot_hyperbola.m
% for hints, see hPositioningTDOACurve.m in LTE System Toolbox example
%
% [input]
% - sensor1_position: [-7, -5]
% - sensor2_position: [5, 1]
% - delta_distance: distance difference between target from sensor 1 and sensor 2
%   positive = target is near sensor 2
%   negative = target is near sensor 1
% [usage]
% pah([-7,-5], [5, 1], 9)

% see https://en.wikipedia.org/wiki/Hyperbola

% set x, y axis ratio to equal, 0 or 1(default, preferable)
axis_equal = 1;

% vector between sensor 1 and sensor 2
delta = sensor1_position - sensor2_position;

% convert cartesian to polar coordinate
[phi, r] = cart2pol(delta(1), delta(2));

% hyperbola parameter
a = -delta_distance/2; % hyperbola vertex, delta distance is equal to distance between two hyperbola vertex
e = r/(2*a);    % eccentricity
c = a * e;
b = sqrt(c^2 - a^2);

% check vertex and focus: focus must be greater than vertex
if abs(a) >= abs(c)
    fprintf('hyperbola vertex(abs(a) = %.2f) must be less than focus(abs(c) = %.2f)\n', abs(a), abs(c));
    return;
end

% hyperbola equation: (x/a)^2 - (y/b)^2 = 1
% hyperbola parametric equation
t = (-3:.01:3)';    % parameter
x = a * cosh(t) - c;
y = b * sinh(t);

[phi2, r2] = cart2pol(x, y);
% [phi2, r2] = cart2pol(real(x), real(y));
[x, y] = pol2cart(phi2 + phi, r2);
x = x + sensor1_position(1);
y = y + sensor1_position(2);

x_max = max(abs(x));
% x_min = min(x);
y_max = max(abs(y));
% y_min = min(y);

figure;
plot(x, y);
if axis_equal
    axis equal;
end
axis([-x_max x_max -y_max y_max]);
grid on;

hold on;
plot(sensor1_position(1), sensor1_position(2), 'ob', ...
    'MarkerSize', 7, 'LineWidth', 2);
plot(sensor2_position(1), sensor2_position(2), 'om', ...
    'MarkerSize', 7, 'LineWidth', 2);
legend('tdoa line', 'sensor 1', 'sensor 2');

end
