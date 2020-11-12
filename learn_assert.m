function [] = learn_assert(ndlrb, nprsrb)
%
% [usage]
% learn_assert(15, 15)
% learn_assert(15, 17)
% 

ndlrb_vec = [6,15,25];

try
    find(ndlrb_vec == ndlrb);
    % ##### below line is NOT needed
    % ##### empty output from "find" command means error
    %     assert(~isempty(idx));
    assert(~(nprsrb > ndlrb));

catch
    disp('wrong ndlrb or nprsrb');
end

end

% prompt = {'enter NDLRB(one of 6,15,25):','enter NPRSRB(1 ~ NDLRB):','enter subframe length(1 or 2)'};
% dlg_title = 'Input';
% num_lines = 1;
% def = {'15','2','1'};
% answer = inputdlg(prompt,dlg_title,num_lines,def);
% S.ndlrb = str2num(answer{1});
% S.nprsrb = str2num(answer{2});
% S.subframe_length = str2num(answer{3});