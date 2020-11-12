function [] = learn_nonlinear_overdetermined

sensor_1 = [4.8541   -0.7114]*1e3;
sensor_2 = [3.9030    4.7709]*1e3;
sensor_3 = [-6.3239    1.8785]*1e3;
sensor_4 = [-4.4540   -7.0001]*1e3;
delta_distance_12 = 1.2491e3;
delta_distance_13 = 1.7176e3;
delta_distance_14 = 3.3571e3;

x0 = [0; 0];  % Make a starting guess at the solution
% x0 = [-5; -5];  % Make a starting guess at the solution
% options = optimoptions('fsolve','Algorithm', 'levenberg-marquardt', 'Display', 'iter', ...
%     'TolFun', 1.0000, 'TolX', 1.0000); % Option to display output
options = optimoptions('fsolve','Algorithm', 'levenberg-marquardt', 'Display', 'iter');
% options = optimoptions('fsolve','Algorithm', 'levenberg-marquardt', 'Display', 'iter');
% options = optimoptions('fsolve','Display','iter'); % Option to display output
f = @(x)tdoa_fun(x, sensor_1, sensor_2, sensor_3, sensor_4, delta_distance_12, delta_distance_13, delta_distance_14);
[x,fval,exitflag] = fsolve(f,x0,options) % Call solver
% [x,fval] = fsolve(@tdoa_fun,x0,options) % Call solver
% [x,fval] = fsolve(@myfun,x0,options) % Call solver


end

function F = myfun(x)
F = [2*x(1) - x(2) - exp(-x(1));
      -x(1) + 2*x(2) - exp(-x(2))];
  
end

function F = tdoa_fun(x, sensor_1, sensor_2, sensor_3, sensor_4, delta_distance_12, delta_distance_13, delta_distance_14)
% sensor_1 = [4.0015e3 -1.3622e3];
% sensor_2 = [2.0857e3 5.5614e3];
% sensor_3 = [-7.7835e3 -2.5652e3];

F = [ ...
    sqrt((sensor_1(1) - x(1))^2 + (sensor_1(2) - x(2))^2) - ...
    sqrt((sensor_2(1) - x(1))^2 + (sensor_2(2) - x(2))^2) + ...
    delta_distance_12;
    
    sqrt((sensor_1(1) - x(1))^2 + (sensor_1(2) - x(2))^2) - ...
    sqrt((sensor_3(1) - x(1))^2 + (sensor_3(2) - x(2))^2) + ...
    delta_distance_13;
    
    sqrt((sensor_1(1) - x(1))^2 + (sensor_1(2) - x(2))^2) - ...
    sqrt((sensor_4(1) - x(1))^2 + (sensor_4(2) - x(2))^2) + ...
    delta_distance_14; ...
    ];

end