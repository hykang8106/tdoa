function [sensor_length, snr_db, randomize_sensor_distance, ...
    use_only_torrieri_method, plot_position, plot_signal] = gui_target_fixed_input_from_uicontrol(S)

lb_sel = get(S.lb(1), {'String', 'value'}); % get sensor length from listbox
% #### tricky!! see comment out part(line 19) in gui_sensor_fixed_input_from_uicontrol.m
sensor_length = str2num(lb_sel{1}{lb_sel{2}});

snr_db = str2num(get(S.ed(1),'String')); % get snr from edit

% get randomize_sensor_distance from checkbox
if get(S.cb(1), 'value')
    randomize_sensor_distance = 0;
else
    randomize_sensor_distance = 1;
end

% get use_only_torrieri_method from checkbox
if get(S.cb(2), 'value')
    use_only_torrieri_method = 1;
else
    use_only_torrieri_method = 0;
end

% get plot_position from checkbox
if get(S.cb(3), 'value')
    plot_position = 1;
else
    plot_position = 0;
end

% get plot_signal from checkbox
if get(S.cb(4), 'value')
    plot_signal = 1;
else
    plot_signal = 0;
end

end