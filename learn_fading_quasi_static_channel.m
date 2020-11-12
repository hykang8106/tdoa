function [] = learn_fading_quasi_static_channel

M = 4;                          % DQPSK modulation order
hMod = comm.DQPSKModulator;       % Create a DPSK modulator
hDemod = comm.DQPSKDemodulator;     % Create a DPSK demodulator

numBits = 10000;                % Each trial uses 10000 bits.
numTrials = 20;                 % Number of BER computations

% Note: In reality, numTrials would be a large number
% to get an accurate estimate of outage probabilities
% or packet error rate.
% Use 20 here just to make the example run more quickly.

% Create Rician channel object.
chan = ricianchan;      % Static Rician channel
chan.KFactor = 3;       % Rician K-factor
% Because chan.ResetBeforeFiltering is 1 by default,
% FILTER resets the channel in each trial below.

% Create an AWGNChannel and ErrorRate calculator System object
hChan = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');
hChan.SNR = 20;
hErrorCalc = comm.ErrorRate;
serVec = zeros(3,numTrials);

% Compute error rate once for each independent trial.
for n = 1:numTrials
    reset(hErrorCalc)
    tx = randi([0 M-1],numBits,1);       % Generate random bit stream
    dpskSig = step(hMod, tx);        % DPSK modulate signal
    fadedSig = filter(chan, dpskSig);    % Apply channel effects
    rxSig = step(hChan,fadedSig);           % Add Gaussian noise.
    rx = step(hDemod,rxSig);       % Demodulate.
    
    % Compute number of symbol errors.
    % Ignore first sample because of DPSK initial condition.
    serVec(:,n) = step(hErrorCalc,tx(2:end),rx(2:end));
end
nErrors = serVec(2,:)
per = mean(nErrors > 0) % Proportion of packets that had errors

end

