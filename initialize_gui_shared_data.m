function [S] = initialize_gui_shared_data

% figure handle
S.fh = zeros(1, 2);
S.fh(:);
length(S.fh);

% text in main figure
S.txmain = 0;

% button group
S.bg = 0;

% radio button
S.rd = zeros(1, 8);
S.rd(:);
length(S.rd);

% push button
S.pb = zeros(1, 2);

% uimenu
S.fm = zeros(1, 11);

% text
S.tx = zeros(1, 8);

% edit
S.ed = zeros(1, 6);

% check box
S.cb = zeros(1, 4);

% list box
S.lb = 0;

S.sensor_position = [];
S.rician_param = [];
S.target_radius_meter = [];
S.uca_radius_meter = [];
S.fading_param_filename = [];
S.rician_param_raw = [];

S.ndlrb = [];
S.nprsrb = [];
S.subframe_length = [];

S.tx_signal = [];
S.fs = [];
S.nfft = [];
S.bw_mhz = [];
S.sample_length = [];

% store current radiobutton
S.SEL = 1;
% store previous radiobutton
S.SEL_old = 1;

% ## not fully tested for small monitor(161128)
root_screen_size = get(0, 'screensize');
% main fig(input) left = 200, main fig(input) width = 850, sub fig(run) width = 350
if root_screen_size(3) > (200 + 850 + 350)
    % for large monitor
    S.second_fig_pos_left = 1070;
else
    % for small monitor 
    S.second_fig_pos_left = 520;
end

end

