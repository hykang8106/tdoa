function [] = learn_hyperbolic_eq_solver(x0)
% solve 2d hyperbolic equation
% ##### only valid when 3 sensor and reference sensor is sensor 1
%
% [input]
% - x0: 2d coordinate, [x y]. initial guess.
% [usage]
% learn_hyperbolic_eq_solver([500 100])

% got from tdoa_position.m
sensor_1 = [4.0015 -1.3622] * 1e3;
sensor_2 = [2.0857 5.5614] * 1e3;
sensor_3 = [-7.7835 -2.5652] * 1e3;
% ### distance difference = (distance between sensor 2 and target) - (distance between sensor 1 and target)
distance_12 = 1.7176e3;
% ### distance difference = (distance between sensor 3 and target) - (distance between sensor 1 and target)
distance_13 = 3.9816e3;

options = optimoptions('fsolve', 'Display', 'iter'); % Option to display output
options
f = @(x)hyperbolic_fun(x, sensor_1, sensor_2, sensor_3, distance_12, distance_13);
[x, fval, exitflag] = fsolve(f, x0, options) % Call solver
% [x,fval] = fsolve(@myfun,x0,options) % Call solver

% #################################
% ### exact solution is [0 0]
% #################################
position_error = sqrt(x(1)^2 + x(2)^2)

end

%%
function F = myfun(x)
F = [2*x(1) - x(2) - exp(-x(1));
      -x(1) + 2*x(2) - exp(-x(2))];
  
end

%%
function F = hyperbolic_fun(x, sensor_1, sensor_2, sensor_3, distance_12, distance_13)

F = [sqrt((sensor_1(1) - x(1))^2 + (sensor_1(2) - x(2))^2) - ...
    sqrt((sensor_2(1) - x(1))^2 + (sensor_2(2) - x(2))^2) + ...
    distance_12;
    
    sqrt((sensor_1(1) - x(1))^2 + (sensor_1(2) - x(2))^2) - ...
    sqrt((sensor_3(1) - x(1))^2 + (sensor_3(2) - x(2))^2) + ...
    distance_13];

end

