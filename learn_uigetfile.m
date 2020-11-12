function [] = learn_uigetfile(filename)
% learn filterspec in uigetfile

if isempty(filename)
    filterspec = 'tdoa_result_sensor_fixed_*.mat';
    [filename, pathname, filterindex] = uigetfile(filterspec);
    if ~filename
        fprintf(2, '######## file selection canceled\n');
        return;
    else
        fprintf('filename = %s\n', filename);
    end
end

% ###############################################################################
% ### if 'filterspec' in uigetfile is set to 'tdoa_result_sensor_fixed_*.mat',
% ### below code is NOT necessary
% ###############################################################################

% I = strfind(filename, 'tdoa_result_sensor_fixed');
% if isempty(I)
%     fprintf(2, '#### filename must have ''tdoa_result_sensor_fixed''\n');
%     return;
% end

end