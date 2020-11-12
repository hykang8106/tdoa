function [] = fm_clear_fading_parameter_call(varargin)

hObject = varargin{3};
S = guidata(hObject);

S.rician_param = [];

uiwait(msgbox({'fading parameter was cleared'},...
    'info','modal'));

guidata(hObject, S);

end

