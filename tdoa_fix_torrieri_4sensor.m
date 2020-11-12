function [] = tdoa_fix_torrieri_4sensor
% (1) matlab nonlinear equation solver: get reference point near target
% (2) torrieri algorithm (see "Statistical Theory of Passive Location Systems" paper written by torrieri)

iteration_length = 5;
% criterion for iteration stop: 
% when absolute difference between initial and estimated is greater than criterion, iteration stop
iteration_criterion = 1; 

c = physconst('LightSpeed');

% ###### set tdoa_position.m input:
% ###### randomize_sensor_distance = 0 (sensor 1 is always reference sensor, emitter is nearest to sensor 1)
% ###### rf_sensor_length = 4
% ###### example: tdoa_position(4, 10, 0)

% ###### got from rf_sensor{i}.Position in tdoa_position.m (line 40)
sensor_1 = [4.8541   -0.7114]*1e3;
sensor_2 = [3.9030    4.7709]*1e3;
sensor_3 = [-6.3239    1.8785]*1e3;
sensor_4 = [-4.4540   -7.0001]*1e3;
% ###### got from delta_distance_meter in tdoa_position.m (line 143)
% ### distance difference = (distance between sensor 2 and target) - (distance between sensor 1 and target)
delta_distance_12 = 1.2491e3;
delta_distance_13 = 1.7176e3;
delta_distance_14 = 3.3571e3;

% x0 = [0; 0];  % Make a starting guess at the solution
% 
% [x_ref, fval, exitflag] = hyperbolic_eq_solver(x0, sensor_1, sensor_2, sensor_3, delta_distance_12, delta_distance_13)
% % eval('[x_ref, fval, exitflag] = hyperbolic_eq_solver(x0, sensor_1, sensor_2, sensor_3, delta_distance_12, delta_distance_13)');
% 
% % ######## considering average of hyperbolic equation solver output #######
% [x_ref, fval, exitflag] = hyperbolic_eq_solver(x0, sensor_1, sensor_3, sensor_4, delta_distance_13, delta_distance_14)

% #### consider average of sensor position
% ############# this option is most reasonable ###########
x_ref = mean([sensor_1; sensor_2; sensor_3; sensor_4])'

% ##### need minus sign #####
Ht = -[delta_distance_12, delta_distance_13, delta_distance_14]' / c;

% sensor dimension: 2 x sensor_length
sensor = [sensor_1' sensor_2' sensor_3' sensor_4'];

% x_ref: column vector (2 x 1)
Rhat = zeros(2, iteration_length);
for n = 1 : iteration_length
    x_est = torrieri(x_ref, sensor, Ht);
    Rhat(:, n) = x_est;
    
    % test iteration criterion
    if abs(sqrt(x_est(1)^2 + x_est(2)^2) - sqrt(x_ref(1)^2 + x_ref(2)^2)) <= iteration_criterion
        break;
    end
    
    x_ref = x_est;
end
Rhat
position_error = sqrt(x_est(1)^2 + x_est(2)^2)

end

%%
function [Rhat] = torrieri(R0, sensor, Ht)
% thanks torrieri
% see "Statistical Theory of Passive Location Systems" paper written by torrieri
%
% [input]
% - R0: initial target position guess, 2 x 1
% - sensor: sensor position, 2 x sensor_length
% - Ht: tdoa between reference sensor and other sensor, (sensor_length - 1) x 1
%   Ht means H * t in eq (81) where t is left term in eq (79) and H is eq (82)
%   in current version, reference sensor is always sensor 1 
% [output]
% - Rhat: estimated target position, 2 x 1

c = 299792458.0;
fs = 3.84e6;
sp = 1 / fs; % sample period
% sp = 1e-7 * sp

bw = 1e6; % signal bandwidth
snr_db = 10;

% see eq (94)
beta_r2 = pi^2 * bw^2 / 3;

% ##### atv: arrival time variance, sigma_t2
% ##### how to compute E and N0? how to related with snr_db?
% see eq (91)
% atv = 1 / ((2 * E / N0) * beta_r2)
% E: energy in received signal, N0/2: two-sided noise power spectral denisty,
% beta_r2: function of bandwidth of signal, see eq (94) for communication signals
% N0 = kT where k is boltzmann constant, T is system noise temperature

% matlab command: systemp, noisepow

% k = physconst('Boltzmann')

% it is nonsense that we set arrival time variance to sample peariod. sample synch must be assumed
% atv = sp;

% (2 * E / N0) = 10
% atv = 1 / (10 * beta_r2);

% ### temporary setting: E/No is assumed to be equal to Eb/No
Eb_over_N0_db = 10;
Eb_over_N0 = 10^(Eb_over_N0_db / 10);
atv = 1 / (2 * Eb_over_N0 * beta_r2);

% get sensor length
sensor_length = size(sensor, 2);

% R0 = repmat(R0, 1, sensor_length);
% repeat R0 column vector to make matrix whose dimension is (2 x sensor_length)
F = (repmat(R0, 1, sensor_length) - sensor)'; % F dimension: sensor_length x 2

D0 = sqrt(F(:, 1).^2 + F(:, 2).^2)
% see eq (78)
F = F ./ repmat(D0, 1, 2)

% for n = 1 : sensor_length
%     F(n, :) = F(n, :) / sqrt(F(n, 1)^2 + F(n, 2)^2);
% end
% F

% ### valid only when sensor 1 is reference sensor
% see eq (82) (modified so that ref sensor is sensor 1)
H = [...
    1 -1 0 0;
    1 0 -1 0;
    1 0 0 -1;
    ];

% atv(sigmat) is related to epsilon(eq (75), (76)), 
% which is arrival time measurement error accounts for
% propagation anomalies, receiver noise, and errors in the assumed station positions.
% see eq (96) (modified so that sensor number is 4)
Ne = [...
    atv 0 0 0;
    0 atv 0 0;
    0 0 atv 0;
    0 0 0 atv;
    ];

% see eq (85)
N = H * Ne * H'

FtHtNi = F' * H' * inv(N)

% see eq (86)
Rhat = R0 + c * inv(FtHtNi * H * F) * FtHtNi * (Ht - H * D0 / c)

% #### P: covariance matrix of estimated target position
% see eq (16),(87)
P = (c^2) * inv(FtHtNi * H * F)

E = eig(P)

% see eq (54)
lambda1 = .5 * [P(1,1) + P(2,2) + sqrt((P(1,1) - P(2,2))^2 + 4 * P(1,2)^2)]
% see eq (55)
lambda2 = .5 * [P(1,1) + P(2,2) - sqrt((P(1,1) - P(2,2))^2 + 4 * P(1,2)^2)]

% GDOP(geometric dilution of precision) is defined as 
% the ratio of the root-mean-square position error Er to the root mean-square ranging error.
% The GDOP indicates how much the fundamental ranging error is magnified 
% by the geometric relation among the transmitter position and the stations.
% see eq (89)
gdop = sqrt(trace(P)) / (c * sqrt(atv))

% CEP(circular error probable) is defined as the radius of the circle 
% that has its center at the mean and contains half the realizations of the random vector.
% CEP is a measure of the uncertainty in the location estimator x relative to its mean E[x].
% see eq (90)
cep = .75 * (c * sqrt(atv)) * gdop
% see eq (74)
cep = .75 * sqrt(P(1,1) + P(2,2))

% coordinate rotation angle in radian, see eq (57)
theta = .5 * atan(2 * P(1,2) / (P(1,1) - P(2,2)))
% coordinate rotation angle in degree, see eq (57)
theta_deg = .5 * atand(2 * P(1,2) / (P(1,1) - P(2,2)))

% concentration ellipse corresponding to probability Pe is defined 
% to be the particular ellipse for which Pe is the probability that x lies inside it. 
% Thus the concentration ellipse is a 2 dimensional measure of accuracy for an unbiased estimator.

% see eq (60)
% k = -2 * log(1 - Pe)

% concentration ellipse: major axes length = 2 * sqrt(k * lambda1), minor axes length = 2 * sqrt(k * lambda2)

% for target position bias, use eq (84), (15)
% see eq (15)
% when linearization error is neglecting, bias is:
% b = inv(FtHtNi * H * F / c^2) * FtHtNi * E[n] 
% where E[n] is expectation of n which is arrival time measurement error
% see eq (75), (76) and see eq (107), (108)

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

