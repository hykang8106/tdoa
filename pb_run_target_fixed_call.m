function [] = pb_run_target_fixed_call(varargin)

hObject = varargin{3};
S = guidata(hObject);

[sensor_length, snr_db, randomize_sensor_distance, ...
    use_only_torrieri_method, plot_position, plot_signal] = ...
    gui_target_fixed_input_from_uicontrol(S);

[position_error_hyperbolic, position_error_linear, position_error_torrieri] = ...
    simulate_tdoa_target_fixed(sensor_length, snr_db, randomize_sensor_distance, ...
    use_only_torrieri_method, plot_position, plot_signal);

guidata(hObject, S);

end

