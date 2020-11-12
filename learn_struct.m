function [] = learn_struct
% learn structure
% #### used to write initialize_gui_shared_data.m

% figure handle, 2 = main_fig_length + sub_fig_length
S.fh = zeros(1, 2);
S.fh(:)
length(S.fh)

S.fh(2) = [];
length(S.fh)
S.fh(:)

% text in main figure
S.txmain = 0;

% button group
S.bg = 0;

% radio button
S.rd = zeros(1, 8);
S.rd(:)
length(S.rd)

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

S.SEL = 1;  % The selected object property of S.bg
S.SEL_old = 1;

% #########################################

s.ed = 10;
s.ed(2) = 20;
s;
length(s);
length(s.ed);

field = 'tx';
value = {1,2};
T = struct(field, value);
T;
length(T);
try
    length(T.tx);
catch
    fprintf(2, '#### error: length(T.tx)\n');
end

field = 'tx';
value = [1,2];
t = struct(field, value);
t;
length(t);
length(t.tx);

field = 'f';
value = {'some text';
         [10, 20, 30];
         magic(5)};
x = struct(field, value);
length(x);
try
    length(x.f);
catch
    fprintf(2, '#### error: length(x.f)\n');
end

y = struct('a',{},'b',{},'c',{});
length(y);

z = struct('a',[],'b',[],'c',[]);
length(z);

end


