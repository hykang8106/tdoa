function [] = learn_fading_compare_result_empirical_theoretical

fs = 10e3;
ts = 1 / fs;
fd = 0;

% Create Rayleigh fading channel object.
chan = rayleighchan(ts, fd);
chan

% Generate data and apply fading channel.
M = 2; % DBPSK modulation order
hMod = comm.DBPSKModulator;       % Create a DPSK modulator
hDemod = comm.DBPSKDemodulator;     % Create a DPSK demodulator
tx = randi([0 M-1],50000,1);        % Generate a random bit stream
dpskSig = step(hMod, tx);       % DPSK modulate the signal
fadedSig = filter(chan,dpskSig);    % Apply the channel effects

% Compute error rate for different values of SNR.
SNR = 0:2:20; % Range of SNR values, in dB.
numSNR = length(SNR);
berVec = zeros(3, numSNR);

% Create an AWGNChannel and ErrorRate calculator System object
hChan = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
hChan
hErrorCalc = comm.ErrorRate;
hErrorCalc

for n = 1:numSNR
    hChan.SNR = SNR(n);
    rxSig = step(hChan,fadedSig);   % Add Gaussian noise
    rx = step(hDemod, rxSig);  % Demodulate
    reset(hErrorCalc)
    % Compute error rate.
    berVec(:,n) = step(hErrorCalc,tx,rx);
end
BER = berVec(1,:);
% Compute theoretical performance results, for comparison.
BERtheory = berfading(SNR,'dpsk',M,1);

% Plot BER results.
figure;
semilogy(SNR,BERtheory,'b-',SNR,BER,'r*');
legend('Theoretical BER','Empirical BER');
xlabel('SNR (dB)'); ylabel('BER');
title('Binary DPSK over Rayleigh Fading Channel');

end
