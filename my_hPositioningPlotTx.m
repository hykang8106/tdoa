%hPositioningPlotTx Plot transmit waveforms
%   hPositioningPlotTx(TX,SR) plots each transmit waveform in the cell
%   array TX at each eNodeB with sampling rate SR.

%   Copyright 2011-2013 The MathWorks, Inc.

function my_hPositioningPlotTx(tx,sr)

    colours = hPositioningColors();
    figure
    hold on;

    t = (0:length(tx{1})-1)/sr;  % Calculate time of samples

    % For each eNodeB plot the transmitted waveform
%     legendstr = cell(1,numel(tx));
    for i = 1:numel(tx)
        plot(t, abs(tx{i}), colours{i}, 'LineWidth', 2);
%         legendstr{i} = sprintf('target');
    end
    grid on;
%     legend(legendstr);
    xlabel('Time (s)');
    ylabel('Absolute value');
    title('Transmitted waveforms(PRS) of target(LTE BS)');

end


