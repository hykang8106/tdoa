function [tx_signal, bw_mhz, nfft] = generate_target_signal_fsk_422mhz(fs, integer_length)

M = 2; 
freq_sep_hz = 7.5e3;
phase_cont = 'discont'; % freq spectrum is more beautiful than 'cont'
% phase_cont = 'cont';

% M = 4; 
% freq_sep_hz = 2.5e3;

% fs constraint: (M - 1) * freq_sep_hz <= fs
if isempty(fs)
    fs = (M - 1) * freq_sep_hz * 2;
end

% determine sample length: sample_length = nsamp * integer_length
% #### if sample length is same as lte prs ndlrb = 15(3840 sample), set integer_length to 480
if isempty(integer_length)
    integer_length = 480;
end

% determine sample length: sample_length = nsamp * integer_length
% #### if sample length is same as lte prs ndlrb = 15(3840 sample), set integer_length to 480
nsamp = 8;

x = randi([0, M - 1], integer_length, 1); % Random signal

tx_signal = fskmod(x, M, freq_sep_hz, nsamp, fs, phase_cont); % Modulate.

bw_mhz = 8.5e-3; % 422mhz signal bandwidth = 8.5khz
% #### dummy, for compatibility with lte prs
nfft = length(tx_signal);

end

