function [] = pb_run_batch_sensor_fixed_call(varargin)

hObject = varargin{3};
S = guidata(hObject);

if isempty(S.sensor_position)
    uiwait(msgbox({'empty sensor position: load sensor postion'},...
        'warning','modal'));
    return;
end

if isempty(S.tx_signal)
    uiwait(msgbox({'empty tx signal: generate tx signal'},...
        'warning','modal'));
    return;
end

% disable all uimenu
set(S.fm,{'enable'},repmat({'off'},length(S.fm),1));

[trial_length, radius_ratio, target_pos_step, snr_db] = ...
    gui_batch_sensor_fixed_input_from_uicontrol(S);

delta_x = target_pos_step;
delta_y = target_pos_step;

S.target_radius_meter = S.uca_radius_meter * radius_ratio;
target_position_limit = S.uca_radius_meter * 3;

x_range = [-1, 1] * S.target_radius_meter;
y_range = [-1, 1] * S.target_radius_meter;

x = x_range(1) : delta_x : x_range(end);
x_length = length(x);
y = y_range(1) : delta_y : y_range(end);
y_length = length(y);

% used to determine whether target is overlapped with sensor
overlap_criterion_meter = 1;

error_torrieri = nan(x_length, y_length);
cep_mean = nan(x_length, y_length);
gdop_mean = nan(x_length, y_length);

plot_signal = 0; plot_position = 0; use_only_torrieri_method = 1;

h = waitbar(0, 'Please wait...', 'name', 'batch program progress', 'windowstyle', 'modal');
steps = x_length * y_length;

u = tic;

for i = 1 : x_length
    for j = 1 : y_length
        
        position_error_torrieri = zeros(1, trial_length);
        % ###### for computing position estimation covariance
        x_est_torrieri = zeros(2, trial_length);
        cep = zeros(1, trial_length);
        gdop = zeros(1, trial_length);
        
        target_position = [x(i), y(j)];
        
        if ~check_target_position_is_good(target_position, S.target_radius_meter)
            fprintf(2, '###### target position is out of range\n');
            return;
        end
        
        if check_target_overlap_sensor(S.sensor_position, target_position, overlap_criterion_meter)
            % target overlap sensor, so skip
            continue;
        end
        
        for k = 1 : trial_length
            [position_error_torrieri(k), x_est_torrieri(:, k), cep(k), gdop(k), ...
                lambda1, lambda2, theta_rad] = ...
                sub_simulate_tdoa_sensor_fixed(S.sensor_position, snr_db, target_position, ...
                plot_signal, plot_position, use_only_torrieri_method, S.target_radius_meter, S.rician_param, ...
                S.tx_signal, S.fs, S.nfft, S.bw_mhz);
        end
        
        idx = position_error_torrieri <= target_position_limit;
        position_error_excl = position_error_torrieri(idx);
        
        mean_torrieri = mean(position_error_excl);
        std_torrieri = std(position_error_excl);
        
        error_torrieri(i, j) = mean_torrieri;
        cep_mean(i, j) = mean(cep(idx));
        gdop_mean(i, j) = mean(gdop(idx));
        
        step = y_length * (i - 1) + j;
        w = toc(u);
        % reference:
        % http://stackoverflow.com/questions/12210583/is-there-a-matlab-function-to-convert-elapsed-seconds-to-hhmmss-format
        z = fix(mod(w, [0, 3600, 60])./[3600, 60, 1]);
        waitbar(step / steps, h, sprintf('%.2f %%, elapsed time = %d : %02d : %02d', ...
            step / steps * 100, z(1), z(2), z(3)));
        
    end
end

close(h);

% for saving result
filename = 'tdoa_result_sensor_fixed.mat';
% #### append date string to filename
[p, name, e] = fileparts(filename);
filename = [name, '_', datestr(now, 'yymmddHHMM'), e];
sensor_position = S.sensor_position;
target_radius_meter = S.target_radius_meter;
target_signal_spec = [S.ndlrb, S.nprsrb, S.subframe_length];
save(filename, 'sensor_position', 'snr_db', 'trial_length', 'x', 'y', 'target_radius_meter', ...
    'delta_x', 'delta_y', 'error_torrieri', 'cep_mean', 'gdop_mean', 'target_signal_spec', ...
    'target_position_limit');

msgbox_str = sprintf('%s was created', filename);
uiwait(msgbox({msgbox_str},'info','modal'));

% enable all uimenu
set(S.fm,{'enable'},repmat({'on'},length(S.fm),1));

guidata(hObject, S);

end

