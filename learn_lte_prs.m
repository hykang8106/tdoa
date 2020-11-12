function [tx_signal] = learn_lte_prs(plot_signal, snr_db, ndlrb, nprsrb, subframe_length)
%
% [input]
% - plot_signal: boolean
% - snr_db: snr in db. if empty, no awgn
% - ndlrb: number of downlink resource block
% - nprsrb: number of position reference signal resource block
% - subframe_length: subframe length. number of prs pulse is proportional to subframe length
% [usage]
% learn_lte_prs(1, [], 15, 10, 1);
% learn_lte_prs(1, 5, 15, 10, 2);

% ndlrb nfft  fs(mhz) sample(per subframe)
% 6	    128   1.92    1920
% 15	256   3.84    3840
% 25	512   7.68    7680
% 50	1024  15.36   15360
% 75	2048  30.72   30720
% 100	2048  30.72   30720

% subcarrier spacing = 15e3 hz
% Choice: 1.4 MHz (6RB) | 3 MHz (15RB) | 5 MHz (25RB) | 10 MHz (50RB) | 15 MHz (75RB) | 20 MHz (100RB)
% http://rfmw.em.keysight.com/wireless/helpfiles/n7624b/Content/Main/Positioning%20Reference%20Signals%20(PRS)%20-%20(Advanced%20LTE%20FDD%20Downlink).htm

% ### [old version] blind touch elephant ###
% ### use get_bw_from_prs_spec_db.m
%
% [ndlrb = 6]
% nprsrb   bw_mhz   pulse_length
% 2        0.4      2
% 4        0.8      2
% 6        1        2
%
% [ndlrb = 15]
% nprsrb   bw_mhz   pulse_length
% 2        0.4      2
% 4        0.8      2
% 6        1        2
% 8        1.5      4
% 10       2        4
% 12       2.4      4
% 14       2.6      4
%
% [ndlrb = 25]
% nprsrb   bw_mhz   pulse_length
% 2        0.4      2
% 4        0.8      2
% 6        1        2
% 8        1.5      4
% 10       2        4
% 12       2.4      4
% 14       2.6      4
% 16       3        4
% 18       3.2      4
% 20       3.4      4
% 22       4        4
% 24       4.2      4

if nprsrb > ndlrb
    error('##### nprsrb MUST NOT be greater than ndlrb');
end

% get lte prs bandwidth from prs specification database
[bw_mhz, fs, nfft, sample_length] = get_bw_from_prs_spec_db(ndlrb, nprsrb, subframe_length);

% generate lte prs for target signal
[tx_signal, fs_enb, nfft_enb] = generate_target_signal_lte_prs(ndlrb, nprsrb, subframe_length);

if fs_enb ~= fs || nfft_enb ~= nfft
    error('##### fs or nfft from prs spec db is not same as enb');
end

if ~isempty(snr_db)
    tx_signal = awgn(tx_signal, snr_db, 'measured', 'db');
end

if plot_signal  
    if ~isempty(snr_db)
        title_text = ...
            sprintf('ndlrb = %d, nprsrb = %d, fs = %f mhz, sub frame = %d, snr = %d db', ...
            ndlrb, nprsrb, fs/1e6, subframe_length, snr_db);
    else
        title_text = ...
            sprintf('ndlrb = %d, nprsrb = %d, fs = %f mhz, sub frame = %d', ...
            ndlrb, nprsrb, fs/1e6, subframe_length);
    end
    
    plot_tx_signal(tx_signal, fs, title_text);
    
    plot_freq_spectrum(tx_signal, fs, title_text);
end

end

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


