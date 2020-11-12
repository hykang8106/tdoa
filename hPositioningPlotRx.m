% hPositioningPlotRx Plot received waveforms
%   hPositioningPlotRx(RX,SR) plots each receive waveform in the cell array
%   RX at each eNodeB with sampling rate SR.

%   Copyright 2011-2013 The MathWorks, Inc.

function hPositioningPlotRx(rx,sr)

    colors = hPositioningColors();
    figure
    hold on;

    t = (0:length(rx{1})-1)/sr;  % Calculate time of samples

    % For each eNodeB plot the received waveform
    legendstr = cell(1, numel(rx));
    for i = 1:numel(rx)
        plot(t, abs(rx{i}), colors{i}, 'LineWidth', 2);
        legendstr{i} = sprintf('enb%d',i);
    end
    legend(legendstr);
    xlabel('Time (s)');
    ylabel('Absolute value');
    title('Received waveforms at UE location');

end
