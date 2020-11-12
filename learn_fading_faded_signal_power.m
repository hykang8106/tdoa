function [] = learn_fading_faded_signal_power

fs = 10e3;
ts = 1 / fs;
fd = 0;

c = rayleighchan(ts, fd);
sig = 1i * ones(2000,1);  % Generate signal
y = filter(c, sig);      % Pass signal through channel
c                       % Display all properties of the channel

% Plot power of faded signal, versus sample number.
figure;
plot(20*log10(abs(y)));

end