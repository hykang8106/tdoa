function [] = fm_generate_lte_prs_signal_call(varargin)

hObject = varargin{3};
S = guidata(hObject);

try
    ndlrb_vec = [6,15,25];
    prompt = {'enter NDLRB(one of 6,15,25):','enter NPRSRB(1 ~ NDLRB):','enter subframe length(1 or 2)'};
    dlg_title = 'Input';
    num_lines = 1;
    def = {'15','2','1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    S.ndlrb = str2num(answer{1});
    S.nprsrb = str2num(answer{2});
    S.subframe_length = str2num(answer{3});
    
    find(ndlrb_vec == S.ndlrb);
    % ##### below line is NOT needed
    % ##### empty output from "find" command means error
    %     assert(~isempty(idx));
    assert(~(S.nprsrb > S.ndlrb));
    
    %             [S.bw_mhz, S.fs, S.nfft, S.sample_length] = get_bw_from_prs_spec_db(S.ndlrb, S.nprsrb, S.subframe_length);
    %
    %             [S.tx_signal, S.fs, S.nfft] = generate_target_signal_lte_prs(S.ndlrb, S.nprsrb, S.subframe_length);
catch
    uiwait(msgbox({'wrong NDLRB or NPRSRB'},'info','modal'));
    %             disp('wrong NDLRB or NPRSRB.');
    return;
end

[S.bw_mhz, S.fs, S.nfft, S.sample_length] = get_bw_from_prs_spec_db(S.ndlrb, S.nprsrb, S.subframe_length);

[S.tx_signal, S.fs, S.nfft] = generate_target_signal_lte_prs(S.ndlrb, S.nprsrb, S.subframe_length);

guidata(hObject, S);

end

