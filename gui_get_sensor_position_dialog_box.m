function [sensor_position, uca_radius_meter, sensor_position_filename] = ...
    gui_get_sensor_position_dialog_box(sensor_position)
% ##### NOT useful: 
% ##### when go to "catch" section, sensor_position is still empty, this cause error

uiwait(msgbox({'empty sensor position', 'dialog box will be opened'},...
    'warning','modal'));
sensor_position_filename = uigetfile('*.txt','select sensor position file');
try
    [sensor_position, uca_radius_meter] = ...
        get_sensor_position_from_file(sensor_position_filename);
catch
    disp('Unable to Load.  Check Name and Try Again.');
end

end

%%
% if isempty(sensor_position)
%     uiwait(msgbox({'empty sensor position', 'dialog box will be opened'},...
%         'warning','modal'));
%     sensor_position_filename = uigetfile('*.txt','select sensor position file');
%     try
%         [sensor_position, uca_radius_meter] = ...
%             get_sensor_position_from_file(sensor_position_filename);
%     catch
%         disp('Unable to Load.  Check Name and Try Again.');
%         return;
%     end
% end