function [] = fm_generate_422mhz_signal_call(varargin)

hObject = varargin{3};
S = guidata(hObject);

% #################### consider to use assert, see learn_assert.m ###################

try
    prompt = {'enter sample rate(Hz):','enter integer length:'};
    dlg_title = 'Input';
    num_lines = 1;
    def = {'0.25e6','480'}; % for default, see learn_fsk_422mhz.m
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    S.fs = str2num(answer{1});
    integer_length = str2num(answer{2});
    
    [S.tx_signal, S.bw_mhz, S.nfft] = generate_target_signal_fsk_422mhz(S.fs, integer_length);
catch
    disp('Unable to Load. Check Name and Try Again.');
end

guidata(hObject, S);

end

