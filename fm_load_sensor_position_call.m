function [] = fm_load_sensor_position_call(varargin)

% ######## DONT USE below line, DONT USE S structure as input of fm_load_sensor_position_call:
% ######## this callback function can destroy uicontrol in "input for sensor fixed" figure
% ######## because all field of S structure is empty(line 1 ~ 61)
%         S = varargin{3};  % Get the structure.

hObject = varargin{3};
S = guidata(hObject);

filterspec = 'sensor_position_*.txt';
S.sensor_position_filename = uigetfile(filterspec,'select sensor position file');
if S.sensor_position_filename == 0
    disp('Unable to Load. Check Name and Try Again.');
    return;
end

[S.sensor_position, S.uca_radius_meter] = ...
    get_sensor_position_from_file(S.sensor_position_filename);

guidata(hObject, S);

%         try
%             [S.sensor_position, S.uca_radius_meter] = ...
%                 get_sensor_position_from_file(S.sensor_position_filename);
%         catch
%             disp('Unable to Load. Check Name and Try Again.')
%         end

end

