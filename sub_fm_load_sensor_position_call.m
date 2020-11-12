function [S] = sub_fm_load_sensor_position_call(S)

filterspec = 'sensor_position_*.txt';
S.sensor_position_filename = uigetfile(filterspec,'select sensor position file');
if S.sensor_position_filename == 0
    disp('Unable to Load. Check Name and Try Again.');
    return;
end

[S.sensor_position, S.uca_radius_meter] = ...
    get_sensor_position_from_file(S.sensor_position_filename);

%         try
%             [S.sensor_position, S.uca_radius_meter] = ...
%                 get_sensor_position_from_file(S.sensor_position_filename);
%         catch
%             disp('Unable to Load. Check Name and Try Again.')
%         end

end
