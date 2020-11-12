function [position_error_hyperbolic, position_error_linear, position_error_torrieri, ...
    x_est_torrieri, cep, gdop, lambda1, lambda2, theta_rad] = ...
    locate_target_tdoa(sensor_position, ref_sensor_number, delta_distance_meter, rx, fs, snr_db, ...
    use_only_torrieri_method, target_position, bw_mhz)
% locate target using tdoa method
%
% ### location method ###
% (1) hyperbolic equation solver
% (2) linear equation solver
% (3) least square estimator(torrieri algorithm)
%
% [input]
% - sensor_position: sensor position. dimension = sensor_length x 2
% - ref_sensor_number: reference sensor number
% - delta_distance_meter: delay distance in meter between reference sensor and other sensor.
%   dimenstion = 1 x sensor_length
% - rx: only for linear_eq_solver function, cell array or matrix having rx signal
%   ##### rx input is NOT nice, but this is temporary for linear equation solver
%   ##### rx input is used for computing time delay using cross correlation of rx signal
% - fs: sample rate
% - snr_db: snr in db
% - use_only_torrieri_method: boolean
% - target_position: target position, dimension must be 2 x 1
% - bw_mhz: target signal(lte prs) bandwidth in mhz

position_error_hyperbolic = [];
position_error_linear = [];
position_error_torrieri = [];
x_est_torrieri = [];
cep = [];
gdop = [];
lambda1 = [];
lambda2 = [];
theta_rad = [];

sensor_length = size(sensor_position, 1);

if use_only_torrieri_method
    % ############ least square estimator(torrieri algorithm)
    [x_est, cep, gdop, lambda1, lambda2, theta_rad] = ...
        thanks_torrieri(sensor_position, ref_sensor_number, delta_distance_meter, snr_db, bw_mhz);
    % x_est_torrieri: dimension = 2 x 1
    x_est_torrieri = x_est;
    position_error_torrieri = sqrt(sum((x_est - target_position).^2));
    fprintf('##### torrieri algorithm: position error = %f meter\n', position_error_torrieri);
else
    if sensor_length == 3
        % ############ least square estimator(torrieri algorithm)
        [x_est, cep, gdop, lambda1, lambda2, theta_rad] = ...
            thanks_torrieri(sensor_position, ref_sensor_number, delta_distance_meter, snr_db, bw_mhz);
        position_error_torrieri = sqrt(sum((x_est - target_position).^2));
        fprintf('##### torrieri algorithm: position error = %f meter\n', position_error_torrieri);
        
        % ############ hyperbolic equation solver
        [x_est, exitflag] = ...
            hyperbolic_equation_solver(sensor_position, ref_sensor_number, delta_distance_meter);
        if exitflag == 1
            x_est
            
            % ######## must be transpose of target_position: target_position'
            position_error_hyperbolic = sqrt(sum((x_est - target_position').^2));
            fprintf('##### hyperbolic equation solver: position error = %f meter\n', position_error_hyperbolic);
        else
            error('##### fsolve error(NOT converged to a root): exit flag = %d\n', exitflag);
            %             fprintf('##### fsolve error(NOT converged to a root): exit flag = %d\n', exitflag);
        end
    elseif sensor_length >= 4
        % ############ least square estimator(torrieri algorithm)
        [x_est, cep, gdop, lambda1, lambda2, theta_rad] = ...
            thanks_torrieri(sensor_position, ref_sensor_number, delta_distance_meter, snr_db, bw_mhz);
        position_error_torrieri = sqrt(sum((x_est - target_position).^2));
        fprintf('##### torrieri algorithm: position error = %f meter\n', position_error_torrieri);
        
        % ############ linear equation solver
        x_est = linear_eq_solver(sensor_position, rx, fs)
        
        position_error_linear = sqrt(sum((x_est - target_position).^2));
        fprintf('##### linear equation solver: position error = %f meter\n', position_error_linear);
    else
        % ### NOT reach here ###
        error('number of sensors must be greater than 2\n');
    end
    
end

end

%%
function [x_est, exitflag] = ...
    hyperbolic_equation_solver(sensor_position, ref_sensor_number, delta_distance_meter)

sensor_length = size(sensor_position, 1);

% ############## DONT delete, use in eval
sensor_1 = sensor_position(1, :);
sensor_2 = sensor_position(2, :);
sensor_3 = sensor_position(3, :);
% initial guess: mean of sensor position
% ############## DONT delete, use in eval
x0 = mean(sensor_position);

other_sensor_number = setdiff(1 : sensor_length, ref_sensor_number);

% ############## DONT delete, use in eval
delta_distance_12 = delta_distance_meter(other_sensor_number(1));
delta_distance_13 = delta_distance_meter(other_sensor_number(2));

command_string = ...
    sprintf('hyperbolic_eq_solver(x0, sensor_%d, sensor_%d, sensor_%d, delta_distance_12, delta_distance_13);', ...
    ref_sensor_number, other_sensor_number(1), other_sensor_number(2));
[x_est, fval, exitflag] = eval(command_string);

end

%% #### DONT delete, use in eval
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

%%
function [T] = linear_eq_solver(sensor_position, rx, fs)
% for tdoa linear equation, see "TDOA_Localization.pdf"(https://github.com/StevenJL/tdoa_localization)
%
% [input]
% - sensor_position: sensor position, dimension = sensor_length x 2
% - rx: cell array having rx signal
% - fs: sample rate
% [output]
% - T: estimated target position. dimension = 2 x 1

v = physconst('LightSpeed');

if iscell(rx)
    % rx was cell array, so convert matrix
    rx = cell2mat(rx);
end
size(rx);
% whos rx

% sort rx signal at each sensor
[Y, I] = sort(sum(abs(rx)), 'descend');
I;
r1 = I(1);
r2 = I(2);
r3 = I(3);

len = size(sensor_position, 1);
timedelayvec = zeros(len, 1);
for i = 1 : len
    timedelayvec(i) = timedelayfunc(rx(:, r1), rx(:, i), fs);
end
timedelayvec;
% ###### when distance of two sensor from target is same and nearest to target,
% ###### delayEst = [1     1    17    38    49]
% ###### this make two element(r1, r2) in timedelayvec become 0
% ###### this make svd error in linear_eq_solver: A(line 177), B(line 179)
% ###### timedelayvec(r2) MUST NOT BE ZERO
if sum(~timedelayvec) == 2
    T = [nan nan]';
    return;
end

A = zeros(len - 2, 1);
B = zeros(len - 2, 1);
C = zeros(len - 2, 1);

x1 = sensor_position(r1, 1);
y1 = sensor_position(r1, 2);

x2 = sensor_position(r2, 1);
y2 = sensor_position(r2, 2);

j = 1;
for i = 1 : len
    if (i == r1) || (i == r2)
        continue;
    end
    
    xi = sensor_position(i, 1);
    yi = sensor_position(i, 2);

    A(j) = (1/(v * timedelayvec(i))) * (-2 * x1 + 2 * xi) - ...
        (1/(v * timedelayvec(r2))) * (-2 * x1 + 2 * x2);
    B(j) = (1/(v * timedelayvec(i))) * (-2 * y1 + 2 * yi) - ...
        (1/(v * timedelayvec(r2))) * (-2 * y1 + 2 * y2);

    Sum1 = (x1^2) + (y1^2) - (xi^2) - (yi^2);
    Sum2 = (x1^2) + (y1^2) - (x2^2) - (y2^2);

    C(j) = v * (timedelayvec(i) - timedelayvec(r2)) + ...
        (1/(v * timedelayvec(i))) * Sum1 - ...
        (1/(v * timedelayvec(r2))) * Sum2;
    
    j = j + 1;
end
A;
B;
C;

M = [A B];

Minv = pinv(M);
T = Minv * (-C);
% x = T(1);
% y = T(2);

end

%%
function out = timedelayfunc(x, y, fs)

c = xcorr(x, y);
[C, I] = max(c);
out = ((length(c) + 1)/2 - I)/fs;
% (length(c) + 1)/2 - I

end

%%
function [x_est, cep, gdop, lambda1, lambda2, theta_rad] = ...
    thanks_torrieri(sensor_position, ref_sensor_number, delta_distance_meter, snr_db, bw_mhz)
%
% [input]
% - sensor_position: sensor_length x 2
% - ref_sensor_number: reference sensor number
% - delta_distance_meter: 1 x sensor_length
% - snr_db: snr in db
% - bw_mhz: bw in mhz
% [output]
% - x_est: dimension = 2 x 1

iteration_length = 5;
% ######## criterion for iteration stop: 
% when absolute difference between initial and estimated is greater than criterion, iteration stop
iteration_criterion = 10;

c = physconst('LightSpeed');

sensor_length = size(sensor_position, 1);

% #### consider average of sensor position
% ############# this option is most reasonable ###########
x_ref = mean(sensor_position)';

% ##### need minus sign #####
Ht = -delta_distance_meter(setdiff(1 : sensor_length, ref_sensor_number))' / c;

[H] = make_torrieri_H_matrix(sensor_length, ref_sensor_number);

% x_ref: column vector (2 x 1)
Rhat = zeros(2, iteration_length);
for n = 1 : iteration_length
    [x_est, lambda1, lambda2, cep, gdop, theta_rad] = ...
        torrieri(x_ref, sensor_position', Ht, H, snr_db, bw_mhz);
    Rhat(:, n) = x_est;
    
    % test iteration criterion
    if abs(sqrt(x_est(1)^2 + x_est(2)^2) - sqrt(x_ref(1)^2 + x_ref(2)^2)) <= iteration_criterion
        break;
    end
    
    x_ref = x_est;
end
Rhat
x_est, lambda1, lambda2, cep, gdop, theta_rad

end

%%
function [Rhat, lambda1, lambda2, cep, gdop, theta_rad] = torrieri(R0, sensor, Ht, H, snr_db, bw_mhz)
% thanks torrieri
% see "Statistical Theory of Passive Location Systems" paper written by torrieri
%
% [input]
% - R0: initial target position guess, 2 x 1
% - sensor: sensor position, 2 x sensor_length
% - Ht: tdoa between reference sensor and other sensor, (sensor_length - 1) x 1
%   Ht means H * t in eq (81) where t is left term in eq (79) and H is eq (82)
%   reference sensor number is not fixed to 1, arbitrary 
% - H:
% - snr_db:
% - bw_mhz:
% [output]
% - Rhat: estimated target position, 2 x 1
% - lambda1: major axis length of concentration ellipse
% - lambda2: minor axis length of concentration ellipse
% - cep: circular error probable
% - gdop: geometric dilution of precision
% - theta_rad: coordinate rotation angle in rad

c = physconst('LightSpeed');
% fs = 3.84e6;
% sp = 1 / fs; % sample period
% sp = 1e-7 * sp

% ###########################################################################################
% ##### LTE PRS(Position Reference Signal) bandwidth: 
% ##### when enb.NDLRB = 15(default), enb.NPRSRB = 2, (see sub_simulate_tdoa_sensor_fixed.m)
% ##### bw = 0.4e6; % 0.4 mhz
% ##### see learn_lte_prs.m
% ###########################################################################################
bw = bw_mhz * 1e6; % signal bandwidth

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

% ########################################################################
% #### more study is NEEDED!!!
% #### temporary setting: E/No is assumed to be equal to Eb/No
% #### temporary setting: snr_db is assumed to be equal to Eb_over_N0_db
% ########################################################################
Eb_over_N0_db = snr_db;
Eb_over_N0 = 10^(Eb_over_N0_db / 10);
atv = 1 / (2 * Eb_over_N0 * beta_r2);

% get sensor length
sensor_length = size(sensor, 2);

% R0 = repmat(R0, 1, sensor_length);
% repeat R0 column vector to make matrix whose dimension is (2 x sensor_length)
F = (repmat(R0, 1, sensor_length) - sensor)'; % F dimension: sensor_length x 2

D0 = sqrt(F(:, 1).^2 + F(:, 2).^2);
% see eq (78)
F = F ./ repmat(D0, 1, 2);

% ### valid only when sensor 1 is reference sensor
% see eq (82) (modified so that ref sensor is sensor 1)
% H = [...
%     1 -1 0 0;
%     1 0 -1 0;
%     1 0 0 -1;
%     ];
% % #### general form for arbitrary sensor length, but reference sensor MUST be sensor 1
% H = diag(-ones(1, sensor_length - 1), 1);
% H = H(1 : end - 1, :); % remove last row
% H(:, 1) = 1; % set first column to 1

% % #### use make_torrieri_H_matrix.m #####
% % ref_sensor_number = 1;
% [H] = make_torrieri_H_matrix(sensor_length, ref_sensor_number);

% atv(sigmat) is related to epsilon(eq (75), (76)), 
% which is arrival time measurement error accounts for
% propagation anomalies, receiver noise, and errors in the assumed station positions.
% see eq (96) (modified so that sensor number is 4)
% Ne = [...
%     atv 0 0 0;
%     0 atv 0 0;
%     0 0 atv 0;
%     0 0 0 atv;
%     ];
Ne = diag(ones(1, sensor_length) * atv); % general form for arbitrary sensor length

% see eq (85)
N = H * Ne * H';

FtHtNi = F' * H' * inv(N);

% see eq (86)
Rhat = R0 + c * inv(FtHtNi * H * F) * FtHtNi * (Ht - H * D0 / c);

% #### P: covariance matrix of estimated target position
% see eq (16),(87)
P = (c^2) * inv(FtHtNi * H * F);

E = eig(P);

% see eq (54)
lambda1 = .5 * [P(1,1) + P(2,2) + sqrt((P(1,1) - P(2,2))^2 + 4 * P(1,2)^2)];
% see eq (55)
lambda2 = .5 * [P(1,1) + P(2,2) - sqrt((P(1,1) - P(2,2))^2 + 4 * P(1,2)^2)];

% GDOP(geometric dilution of precision) is defined as 
% the ratio of the root-mean-square position error Er to the root mean-square ranging error.
% The GDOP indicates how much the fundamental ranging error is magnified 
% by the geometric relation among the transmitter position and the stations.
% see eq (89)
gdop = sqrt(trace(P)) / (c * sqrt(atv));

% CEP(circular error probable) is defined as the radius of the circle 
% that has its center at the mean and contains half the realizations of the random vector.
% CEP is a measure of the uncertainty in the location estimator x relative to its mean E[x].
% see eq (74)
cep = .75 * sqrt(P(1,1) + P(2,2));
% see eq (90)
cep = .75 * (c * sqrt(atv)) * gdop;

% coordinate rotation angle in radian, see eq (57)
theta_rad = .5 * atan(2 * P(1,2) / (P(1,1) - P(2,2)));
% coordinate rotation angle in degree, see eq (57)
theta_deg = .5 * atand(2 * P(1,2) / (P(1,1) - P(2,2)));

% concentration ellipse corresponding to probability Pe is defined 
% to be the particular ellipse for which Pe is the probability that x lies inside it. 
% Thus the concentration ellipse is a 2 dimensional measure of accuracy for an unbiased estimator.

% concentration ellipse corresponding to probability Pe: 
% major axes length = 2 * sqrt(k * lambda1) 
% minor axes length = 2 * sqrt(k * lambda2)
% where k = -2 * log(1 - Pe)
% see eq (60)

% for target position bias, use eq (84), (15)
% see eq (15)
% when linearization error is neglecting, bias is:
% b = inv(FtHtNi * H * F / c^2) * FtHtNi * E[n] 
% where E[n] is expectation of n which is arrival time measurement error
% see eq (75), (76) and see eq (107), (108)

% estimator target position(Rhat) accuracy:
% rms error, epsilon_r
% epsilon_r^2 = trace(P) + sum(bias.^2)
% where P is covariance matrix of estimated target position
% see eq (51), (52)

end

