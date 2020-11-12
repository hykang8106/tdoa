function [] = learn_lte_prs(plot_signal, ndlrb, nprsrb, snr_db)
%
% [input]
% - plot_signal:
% - ndlrb: number of downlink resource block
% - nprsrb: number of position reference signal resource block
% - snr_db:
% [usage]
% learn_lte_prs(1, 15, 10, [])

% ndlrb nfft  fs(mhz) sample
% 6	    128   1.92    1920
% 15	256   3.84    3840
% 25	512   7.68    7680
% 50	1024
% 75	2048
% 100	2048

% ### blind touch elephant ###
%
% [ndlrb = 6]
% nprsrb   bw_mhz   pulse_number
% 2        0.4      2
% 4        0.8      2
% 6        1        2
%
% [ndlrb = 15]
% nprsrb   bw_mhz   pulse_number
% 2        0.4      2
% 4        0.8      2
% 6        1        2
% 8        1.5      4
% 10       2        4
% 12       2.4      4
% 14       2.6      4
%
% [ndlrb = 25]
% nprsrb   bw_mhz   pulse_number
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

% target emitter to be located is lte base station
enb = lteRMCDL('R.5');   % Get configuration based on RMC
enb.NCellID = 0;     % Set unique cell identity
enb.TotSubframes = 1;  % Number of subframes to generate
enb.NPRSRB = nprsrb;        % Bandwidth of PRS in resource blocks
% enb.NPRSRB = 2;        % Bandwidth of PRS in resource blocks, default = 2
enb.IPRS = 0;          % PRS configuration index
enb.PRSPeriod = 'On';  % PRS present in all subframes
enb.NDLRB = ndlrb;  % default = 15
% enb.CyclicPrefix = 'Extended';
% enb.DuplexMode = 'TDD';
% Warning: Using default value for parameter field TDDConfig (0) 
% Warning: Using default value for parameter field SSC (0) 
enb

% get ofdm modulation related information
info = lteOFDMInfo(enb);
info

fs = info.SamplingRate;

resource_grid = lteDLResourceGrid(enb);        % Empty resource grid
% resource_grid = repmat(resource_grid, 1, 2);
x = ltePRS(enb) % Map PRS REs
y = ltePRSIndices(enb)
resource_grid(y) = x;
% resource_grid(ltePRSIndices(enb)) = ltePRS(enb); % Map PRS REs
tx_signal = lteOFDMModulate(enb, resource_grid);        % OFDM modulate
size(resource_grid)
tx_signal_length = size(tx_signal) % column vector

if ~isempty(snr_db)
    tx_signal = awgn(tx_signal, snr_db, 'measured', 'db');
    title_text = sprintf('ndlrb = %d, nprsrb = %d, fs = %f mhz, snr = %d db', ndlrb, nprsrb, fs/1e6, snr_db);
else
    title_text = sprintf('ndlrb = %d, nprsrb = %d, fs = %f mhz', ndlrb, nprsrb, fs/1e6);
end

% Plot transmit waveforms from target emitter
if plot_signal
    plot_tx_signal(tx_signal, fs);
    
    plot_freq_spectrum(tx_signal, fs, title_text);
end

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
title(title_text);

% Hs = spectrum.periodogram;
% % Hs = spectrum.welch;
% psd(Hs, signal, 'Fs', fs, 'CenterDC', true);
% % psd(Hs, signal, 'Fs', fs, 'SpectrumType','twosided');

%     % create periodogram object for spectrum plot
%     hp = spectrum.periodogram('Hamming');
%     
%     % plot spectrum of signal
%     subplot(2,1,1);
%     psd_hp = psd(hp, x, 'Fs', fs, 'NFFT', length(x), 'CenterDC', true);
%     plot(psd_hp);
%     ylabel('PSD in dB/Hz');
%     first_text = sprintf('[spectrum] %s', title_text);
%     title(first_text, 'interpreter', 'none');

end


