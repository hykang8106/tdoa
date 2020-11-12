function [trial_length, radius_ratio, target_pos_y, target_pos_x, ...
    snr_db, subframe_length, nprsrb, ndlrb] = gui_sensor_fixed_input_from_uicontrol(S)

trial_length = str2num(get(S.ed(1),'String')); % trial length
% get 'string', 'value' property of edit uicontrol
get(S.ed(1), {'string', 'value'});
% get output example = '100' [0], ##### what is [0]? [ans] no meaning in edit uicontrol, index begin from 1

radius_ratio = str2num(get(S.ed(2),'String')); % radius ratio

target_pos_y = str2num(get(S.ed(3),'String')); % y of target position

target_pos_x = str2num(get(S.ed(4),'String')); % x of target position

snr_db = str2num(get(S.ed(5),'String')); % snr in db

% ca = cell array
ca = get(S.pp(1), {'string', 'value'}); % subframe length
subframe_length = str2num(ca{1}{ca{2}});

ca = get(S.pp(2), {'string', 'value'}); % nprsrb
nprsrb = str2num(ca{1}{ca{2}});

ca = get(S.pp(3), {'string', 'value'}); % ndlrb
ndlrb = str2num(ca{1}{ca{2}});

end