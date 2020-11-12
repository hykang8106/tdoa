function [] = tdoa_fix_torrieri
% (1) matlab nonlinear equation solver: get reference point near target
% (2) torrieri algorithm (see "Statistical Theory of Passive Location Systems" paper written by torrieri)

c = 299792458.0;

% set tdoa_position.m input:
% randomize_sensor_distance = 0 (sensor 1 is always reference sensor, emitter is nearest to sensor 1)
% rf_sensor_length = 4
sensor_1 = [4.8541   -0.7114]*1e3;
sensor_2 = [3.9030    4.7709]*1e3;
sensor_3 = [-6.3239    1.8785]*1e3;
sensor_4 = [-4.4540   -7.0001]*1e3;

delta_distance_12 = 1.2491e3;
delta_distance_13 = 1.7176e3;
delta_distance_14 = 3.3571e3;

x0 = [0; 0];  % Make a starting guess at the solution

[x_ref, fval, exitflag] = hyperbolic_eq_solver(x0, sensor_1, sensor_2, sensor_3, delta_distance_12, delta_distance_13)
% eval('[x_ref, fval, exitflag] = hyperbolic_eq_solver(x0, sensor_1, sensor_2, sensor_3, delta_distance_12, delta_distance_13)');

% % ######## considering average of hyperbolic equation solver output #######
% [x_ref, fval, exitflag] = hyperbolic_eq_solver(x0, sensor_1, sensor_3, sensor_4, delta_distance_13, delta_distance_14)

% #### consider average of sensor position
x_ref = mean([sensor_1; sensor_2; sensor_3; sensor_4])'

% ##### need minus sign #####
Ht = -[delta_distance_12, delta_distance_13, delta_distance_14]' / c;

% sensor dimension: 2 x sensor_length
sensor = [sensor_1' sensor_2' sensor_3' sensor_4'];

% x_ref: column vector (2 x 1)
Rhat = torrieri(x_ref, sensor, Ht)

end

%%
function [Rhat] = torrieri(R0, sensor, Ht)
% thanks torrieri
% see "Statistical Theory of Passive Location Systems" paper written by torrieri
%
% [input]
% - R0: 2 x 1
% - sensor: 2 x sensor_length
% - Ht: (sensor_length - 1) x 1
% [output]
% - Rhat: 2 x 1

c = 299792458.0;
fs = 3.84e6;
sp = 1 / fs; % sample period
% sp = 1e-7 * sp

bw = 1e6; % signal bandwidth
snr_db = 10;

% eq (94)
beta_r2 = pi^2 * bw^2 / 3;

% ##### atv: arrival time variance
% ##### how to compute E and N0? how to related with snr_db?
% eq (91)
% atv = 1 / ((2 * E / N0) * beta_r2)

% it is nonsense that we set arrival time variance to sample peariod. sample synch must be assumed
% atv = sp;

% (2 * E / N0) = 10
atv = 1 / (10 * beta_r2);

% get sensor length
sensor_length = size(sensor, 2);

% R0 = repmat(R0, 1, sensor_length);
% repeat R0 column vector to make matrix whose dimension is (2 x sensor_length)
F = (repmat(R0, 1, sensor_length) - sensor)'; % F dimension: sensor_length x 2

D0 = sqrt(F(:, 1).^2 + F(:, 2).^2)
% eq (78)
F = F ./ repmat(D0, 1, 2)

% for n = 1 : sensor_length
%     F(n, :) = F(n, :) / sqrt(F(n, 1)^2 + F(n, 2)^2);
% end
% F

% ### valid only when sensor 1 is reference sensor
% eq (82) (modified so that ref sensor is sensor 1)
H = [...
    1 -1 0 0;
    1 0 -1 0;
    1 0 0 -1;
    ];

% eq (96) (modified so that sensor number is 4)
Ne = [...
    atv 0 0 0;
    0 atv 0 0;
    0 0 atv 0;
    0 0 0 atv;
    ];

% eq (85)
N = H * Ne * H'

FtHtNi = F' * H' * inv(N)

% eq (86)
Rhat = R0 + c * inv(FtHtNi * H * F) * FtHtNi * (Ht - H * D0 / c)

% eq (87)
P = (c^2) * inv(FtHtNi * H * F)

E = eig(P)

% eq (54)
lambda1 = .5 * [P(1,1) + P(2,2) + sqrt((P(1,1) - P(2,2))^2 + 4 * P(1,2)^2)]
% eq (55)
lambda2 = .5 * [P(1,1) + P(2,2) - sqrt((P(1,1) - P(2,2))^2 + 4 * P(1,2)^2)]

% eq (89)
gdop = sqrt(trace(P)) / (c * sqrt(atv))

% eq (90)
cep = .75 * (c * sqrt(atv)) * gdop
% eq (74)
cep = .75 * sqrt(P(1,1) + P(2,2))

% coordinate rotation angle, eq (57)
theta = .5 * atan(2 * P(1,2) / (P(1,1) - P(2,2)))

theta_deg = .5 * atand(2 * P(1,2) / (P(1,1) - P(2,2)))

end

%%
function [x, fval, exitflag] = ...
    hyperbolic_eq_solver(x0, sensor_1, sensor_2, sensor_3, delta_distance_12, delta_distance_13)

options = optimoptions('fsolve','Display','off'); % Option to display output
f = @(x)hyperbolic_fun(x, sensor_1, sensor_2, sensor_3, delta_distance_12, delta_distance_13);
[x, fval, exitflag] = fsolve(f,x0,options); % Call solver

end

%%
function F = hyperbolic_fun(x, sensor_1, sensor_2, sensor_3, delta_distance_12, delta_distance_13)

F = [sqrt((sensor_1(1) - x(1))^2 + (sensor_1(2) - x(2))^2) - ...
    sqrt((sensor_2(1) - x(1))^2 + (sensor_2(2) - x(2))^2) + ...
    delta_distance_12;
    sqrt((sensor_1(1) - x(1))^2 + (sensor_1(2) - x(2))^2) - ...
    sqrt((sensor_3(1) - x(1))^2 + (sensor_3(2) - x(2))^2) + ...
    delta_distance_13];

end

