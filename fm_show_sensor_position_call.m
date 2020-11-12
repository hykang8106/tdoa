function [] = fm_show_sensor_position_call(varargin)

hObject = varargin{3};
S = guidata(hObject);

if isempty(S.sensor_position)
    uiwait(msgbox({'empty sensor position: load sensor postion'},...
        'warning','modal'));
    return;
end

title_text = 'sensor position';
% ####### radius_ratio: temporary #######
% ####### radius_ratio determine in pb_sensor_fixed_call, so here it have no value
radius_ratio = 1.5;
S.target_radius_meter = S.uca_radius_meter * radius_ratio;
plot_sensor_position_only(S.sensor_position, S.target_radius_meter, title_text, ...
    S.fading_param_filename);

guidata(hObject, S);

end

