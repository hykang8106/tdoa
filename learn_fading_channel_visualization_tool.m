function [] = learn_fading_channel_visualization_tool

% Three-Path Rayleigh channel
h = rayleighchan(1/100000, 130, [0 1.5e-5 3.2e-5], [0, -3, -3]);
tx = randi([0 1],500,1);          % Random bit stream
hmod = comm.DBPSKModulator;       % Create DBPSKModulator
dpskSig = step(hmod,tx);          % DPSK signal
h.StoreHistory = true;            % Allow states to be stored
y = filter(h, dpskSig);           % Run signal through channel
plot(h);                          % Call Channel Visualization Tool

end
