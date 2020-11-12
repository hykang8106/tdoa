function [] = fm_show_fading_parameter_call(varargin)

hObject = varargin{3};
S = guidata(hObject);

if isempty(S.rician_param)
    uiwait(msgbox({'fading parameter is empty'},...
        'info','modal'));
    return;
end

sensor_length = size(S.rician_param_raw, 1);
rnames = cell(1, sensor_length);
for n = 1 : sensor_length
    rnames{n} = sprintf('sensor%d', n);
end

f = figure('Position',[200 200 450 200], ...
    'menubar','none',...
    'name','fading parameters',...
    'numbertitle','off',...
    'resize','off');
dat = S.rician_param_raw;
cnames = {'K factor','path delay(delta ratio)','path loss(dB)'};
t = uitable('Parent',f,'Data',dat,'ColumnName',cnames,...
    'RowName',rnames,'Position',[20 20 370 170]);

guidata(hObject, S);

end

