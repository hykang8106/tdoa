function [] = learn_fading_channel_delay

M = 2;                          % DQPSK modulation order
bitRate = 50000;
hMod = comm.DBPSKModulator;       % Create a DPSK modulator
hDemod = comm.DBPSKDemodulator;     % Create a DPSK demodulator

% Create Rayleigh fading channel object.
ch = rayleighchan(1/bitRate,4,[0 0.5/bitRate],[0 -10]);
delay = ch.ChannelFilterDelay;

tx = randi([0 M-1],50000,1);        % Generate random bit stream
dpskSig = step(hMod,tx);        % DPSK modulate signal
fadedSig = filter(ch,dpskSig);      % Apply channel effects
rx = step(hDemod,fadedSig);   % Demodulate signal

% Compute bit error rate, taking delay into account.
hErrorCalc = comm.ErrorRate('ReceiveDelay', delay);
berVec = step(hErrorCalc,tx,rx);
ber = berVec(1)
num = berVec(2)

end
