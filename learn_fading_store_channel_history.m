function [] = learn_fading_store_channel_history

h = rayleighchan(1/100000, 130);  % Rayleigh channel
tx = randi([0 1],10,1);           % Random bit stream
hmod = comm.DBPSKModulator;       % Create DBPSK Modulator
dpskSig = step(hmod,tx);          % Process data by calling the step method
h.StoreHistory = true;            % Allow states to be stored
y = filter(h, dpskSig);           % Run signal through channel
h.PathGains                       % Display the stored path gains data

end