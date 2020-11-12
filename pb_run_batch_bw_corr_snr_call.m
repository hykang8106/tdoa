function [] = pb_run_batch_bw_corr_snr_call(varargin)

hObject = varargin{3};
S = guidata(hObject);

if isempty(S.sensor_position)
    uiwait(msgbox({'empty sensor position: load sensor postion'},...
        'warning','modal'));
    return;
end

target_position = str2num(char(get(S.ed(1),'String')));
radius_ratio = str2num(get(S.ed(2),'String'));
trial_length = str2num(get(S.ed(3),'String'));
snr_vec = str2num(char(get(S.ed(4),'String')));
nprsrb_vec = str2num(char(get(S.ed(5),'String')));
nsubframe_vec = str2num(char(get(S.ed(6),'String')));

lb_sel = get(S.lb(1), {'String', 'value'}); % get sensor length from listbox
% #### tricky!! see comment out part(line 19) in gui_sensor_fixed_input_from_uicontrol.m
ndlrb = str2num(lb_sel{1}{lb_sel{2}});

if nprsrb_vec(end) > ndlrb
    uiwait(msgbox({'NPRSRB(end) MUST NOT be greater than NDLRB'},...
        'warning','modal'));
    return;
    %                 fprintf(2, '##### nprsrb MUST NOT be greater than ndlrb\n');
    %                 return;
end

[filename] = sub_batch_simulate_tdoa_bw_corr_snr(S.sensor_position, S.uca_radius_meter, ...
    target_position, radius_ratio, trial_length, S.rician_param, ...
    ndlrb, snr_vec, nprsrb_vec, nsubframe_vec);

msgbox_str = sprintf('%s was created', filename);
uiwait(msgbox({msgbox_str},'info','modal'));

guidata(hObject, S);

end

