%hPositioningPlotCorr Plot eNodeB correlation
%   hPositioningPlotCorr(CORR,SR) plots the correlation of each eNodeB
%   in CORR given a sampling rate SR.

%   Copyright 2011-2013 The MathWorks, Inc.

function hPositioningPlotCorr(corr,sr)

    colors = hPositioningColors();
    figure;
    hold on;
    
    % Plot correlation for each eNodeB
    t=(0:length(corr{1})-1)/sr;
    legendstr = cell(1, numel(corr));
    for i=1:numel(corr)
        plot(t, abs(corr{i}), colors{i}, 'LineWidth', 2);
        legendstr{i} = sprintf('enb%d', i);
    end
    legend(legendstr);
    xlabel('Time (s)');
    ylabel('Absolute value');
    title('Receiver correlations at UE location');
    
    % Plot markers are correlation peaks
    for i=1:numel(corr)
        c=abs(corr{i});
        j=find(c==max(c), 1);
        plot(t(j), c(j),strcat(colors{i}, 'o'), 'MarkerSize', 5, ...
            'MarkerFaceColor', colors{i}, 'LineWidth', 2);
        legendstr{i} = sprintf('enb%d', i);
    end

end
