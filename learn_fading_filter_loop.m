function [] = learn_fading_filter_loop

% Set up parameters.
M = 4;              % QPSK modulation order
hMod = comm.QPSKModulator;
bitRate = 50000;    % Data rate is 50 kb/s
numTrials = 1250;    % Number of iterations of loop

% Create Rayleigh fading channel object.
ch = rayleighchan(1/bitRate,4,[0 2e-5],[0 -9]);
% Indicate that FILTER should not reset the channel
% in each iteration below.
ch.ResetBeforeFiltering = 0;

% Initialize scatter plot.
scatterPlot = commscope.ScatterPlot;

% Apply channel in a loop, maintaining continuity.
% Plot only the current data in each iteration.
for n = 1:numTrials
   tx = randi([0 M-1],500,1);     % Generate random bit stream
   pskSig = step(hMod,tx);         % PSK modulate signal
   fadedSig = filter(ch, pskSig); % Apply channel effects

   % Plot the new data from this iteration.
   update(scatterPlot,fadedSig);
end

end