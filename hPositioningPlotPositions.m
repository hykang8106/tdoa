%hPositioningPlotPositions Plot eNodeB positions
%   hPositioningPlotPositions(ENB) plots the positions of eNodeBs in the
%   cell array ENB.

%   Copyright 2011-2013 The MathWorks, Inc.

function hPositioningPlotPositions(enb)

    colours = hPositioningColors();
    styles = hStyles();
    
    figure;
    hold on;
    legendstr = cell(1, numel(enb));
    
    % Plot the position of each eNodeB
    for i = 1:numel(enb)
        plot(enb{i}.Position(1),enb{i}.Position(2), ...
            strcat(colours{i}, styles{i}), ...
            'MarkerSize', 7, 'LineWidth', 2);
        legendstr{i} = sprintf('enb%d',i);
    end
    
    % Plot an origin point
    plot(0, 0, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'c', ...
        'LineWidth', 2);
    
    axis([-11000 11000 -11000 11000]);
    legend([legendstr 'UE']);
    xlabel('X position (m)');
    ylabel('Y position (m)');
    title('Positions');

end

%hStyles Plotting styles
%  STYLES = hStyles() returns a lookup table of styles for plotting.

%   Copyright 2011-2013 The MathWorks, Inc.

function styles = hStyles()

    styles = {'*', 'd', 's', 'p', '^'};

end
