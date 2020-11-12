function [] = learn_fsk_422mhz(snr_db, plot_signal, fs, integer_length)
% learn fsk to generate 422mhz target signal whose bandwidth is 8.5khz
% #### used in tdoa simulation for 422mhz radio signal(narrow bandwidth = 8.5 khz)
% #### power of 422mhz radio station may be 5W, tx_signal amplitude is right for 5W?
%
% [input]
% - snr_db:
% - plot_signal:
% - fs: sample rate. fs constraint: (M-1)*freq_sep_hz <= fs
%   if empty, set to 15e3 hz
% - integer_length: integer length, determine sample length
%   if empty, set to 480
%   if sample length is same as lte prs ndlrb = 15(3840 sample), set integer_length to 480
%
% [usage]
% learn_fsk_422mhz(5, 1, [], [])
% learn_fsk_422mhz(5, 1, 15e3, 480)
% learn_fsk_422mhz(5, 1, 16e3, 1000)
% 

% ####### tx_signal length = 8 * integer_length
% ####### integer_length >= 240
% ####### fs >= 15e3
% ####### bw_mhz always return with 8.5e-3(422mhz signal bandwidth = 8.5khz)
% ####### nfft: length(tx_signal), dummy, for compatibility with lte prs
[tx_signal, bw_mhz, nfft] = generate_target_signal_fsk_422mhz(fs, integer_length);
sample_length = length(tx_signal);

if ~isempty(snr_db)
    tx_signal = awgn(tx_signal, snr_db, 'measured', 'db');
end

if plot_signal
    
    title_text = sprintf('snr = %d db, fs = %f khz, sample = %d', snr_db, fs / 1e3, sample_length);

    plot_tx_signal(tx_signal, fs, title_text);
    
    plot_freq_spectrum(tx_signal, fs, title_text);
    
end

end

%%
% function [tx_signal, bw_mhz, nfft] = generate_fsk_422mhz(fs, integer_length)
% 
% M = 2; 
% freq_sep_hz = 7.5e3;
% phase_cont = 'discont'; % freq spectrum is more beautiful than 'cont'
% % phase_cont = 'cont';
% 
% % M = 4; 
% % freq_sep_hz = 2.5e3;
% 
% % fs constraint: (M - 1) * freq_sep_hz <= fs
% if isempty(fs)
%     fs = (M - 1) * freq_sep_hz * 2;
% end
% 
% % determine sample length: sample_length = nsamp * integer_length
% % #### if sample length is same as lte prs ndlrb = 15(3840 sample), set integer_length to 480
% if isempty(integer_length)
%     integer_length = 480;
% end
% 
% % determine sample length: sample_length = nsamp * integer_length
% % #### if sample length is same as lte prs ndlrb = 15(3840 sample), set integer_length to 480
% nsamp = 8;
% 
% x = randi([0, M - 1], integer_length, 1); % Random signal
% 
% tx_signal = fskmod(x, M, freq_sep_hz, nsamp, fs, phase_cont); % Modulate.
% 
% bw_mhz = 8.5e-3; % 422mhz signal bandwidth = 8.5khz
% % #### dummy, for compatibility with lte prs
% nfft = length(tx_signal);
% 
% end

%%
function plot_tx_signal(tx_signal, fs, title_text)
% ############## MUST BE LOCAL FUNCTION
% ############## matlab file whose name is same as this function exist in tdoa directory

figure;

t = (0 : length(tx_signal) - 1) / fs;  % Calculate time of samples
plot(t, abs(tx_signal));
grid on;
xlabel('time in sec');
ylabel('absolute value');
title(['[time domain] ', title_text]);

end

%%
function [] = plot_freq_spectrum(signal, fs, title_text)

sig_len = length(signal);

window = hamming(sig_len);
nfft = sig_len;
freqrange = 'centered';
spectrumtype = 'power';
% spectrumtype = 'psd';
[pxx, f] = periodogram(signal, window, nfft, fs, freqrange, spectrumtype);

figure;
plot(f / 1e6, 10*log10(pxx));
xlabel('freq in mhz');
ylabel('power in dbm');
grid on;
title(['[freq domain] ', title_text]);

end

