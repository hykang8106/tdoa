function [] = fm_load_fading_parameter_call(varargin)

hObject = varargin{3};
S = guidata(hObject);

% ###### uigetfile move to "try" section is better?
S.fading_param_filename = uigetfile('*.xlsx','select fading parameter file');
if S.fading_param_filename == 0
    disp('Unable to Load. Check Name and Try Again.');
    return;
end

try
    prompt = {'Enter Excel Range:'};
    dlg_title = 'Input';
    num_lines = 1;
    def = {'b3:d7'};
    xlsrange = inputdlg(prompt,dlg_title,num_lines,def); % inputdlg output is 'cell'
    xlsrange = cell2mat(xlsrange); % xlsrange input of xlsread MUST BE 'char', not 'cell'
    [S.rician_param, S.rician_param_raw] = ...
        get_rician_parameter_from_excel_file(S.fading_param_filename, xlsrange);
    S.rician_param;
    %             if size(S.rician_param, 1) ~= size(S.sensor_position, 1)
    %                 error('##### row length of rician parameter must be same as sensor length\n');
    %             end
catch
    disp('Unable to Load. Check Name and Try Again.');
    return;
end

if size(S.rician_param, 1) ~= size(S.sensor_position, 1)
    errordlg('##### row length of rician parameter must be same as sensor length', ...
        'sensor length error', 'modal');
    return;
end

guidata(hObject, S);

end

