function [faded_signal, chan] = rician_fading(tx_signal, fs, k, tau, pdb)

% ###################################################################################################
% ### Setting the maximum path delay greater than 100 samples may generate an ¡®Out of memory' error.
% ### when fs = 3.84e6 hz(sample period = 0.26 us), max path delay MUST BE LESS than 26 us
% ###################################################################################################

ts = 1 / fs;
fd = 0; % max doppler shift
% ##### rician k: ratio between specular power and diffuse power for a direct los path, 
% k = 1 ~ 10, rayleigh fading: k = 0, default: k = 1
% k = 1;
% k = 1;
% k = [3, 1];

% ########## For outdoor environments, 
% ########## path delays after the first are typically between 100 ns(= 1e-7 s) and 10 us(= 1e-5 s)
% ########## 100 ns = 30 m, 10 us = 3 km
% tau = [0, 7e-6];
% pdb = [-20, -30];
% pdb = [0, -3];
chan = ricianchan(ts, fd, k, tau, pdb);
chan;
% chan.PathDelays = 0;
% % chan.PathDelays = delay;
% chan.AvgPathGaindB = 0;
% chan.AvgPathGaindB = -20;
% chan.AvgPathGaindB = -PLdB;
% chan.NormalizePathGains = 0;
% ###############################################################################
% ##### StoreHistory must be false if MaxDopplerShift is zero.
% ##### channel visualization tool CANNOT BE USED when MaxDopplerShift is zero.
% ###############################################################################
% chan.StoreHistory = true; 
% chan

% % Ts = 1e-4;
% % fd = 100;
% k = [3 1];
% tau = [0 1e-5 1.5e-5 3e-5];
% pdb = [0 -1 -2 -2.5];
% h = ricianchan(Ts, fd, k, tau, pdb);

faded_signal = filter(chan, tx_signal);
chan;
size(faded_signal);

if chan.ChannelFilterDelay
    faded_signal = [faded_signal(chan.ChannelFilterDelay + 1 : end); zeros(chan.ChannelFilterDelay, 1)];
end

end

