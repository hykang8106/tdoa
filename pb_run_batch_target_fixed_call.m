function [] = pb_run_batch_target_fixed_call(varargin)

hObject = varargin{3};
S = guidata(hObject);

% ####################################################################
% #### MUST USE S.ed(1)
% #### DONT USE S.ed because S.ed is edit array!!
% ####################################################################

% #### str2num input MUST BE 'char array', NOT cell
% #### so char function is used to convert 'cell' to 'char array'
sensor_length = str2num(char(get(S.ed(1),'String')));
snr_db = str2num(char(get(S.ed(2),'String')));
trial_length = str2num(get(S.ed(3),'String'));

filename = batch_simulate_tdoa_target_fixed(sensor_length, snr_db, trial_length);

msgbox_str = sprintf('%s was created', filename);
uiwait(msgbox({msgbox_str},'info','modal'));

guidata(hObject, S);

end

