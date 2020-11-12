function plot_tx_signal(tx_signal, fs)

colours = hPositioningColors();
figure
hold on;

t = (0:length(tx_signal)-1)/fs;  % Calculate time of samples

plot(t, abs(tx_signal), colours{1}, 'LineWidth', 2);

% For each eNodeB plot the transmitted waveform
%     legendstr = cell(1,numel(tx));
% for i = 1:numel(tx)
%     plot(t, abs(tx{i}), colours{i}, 'LineWidth', 2);
%     %         legendstr{i} = sprintf('target');
% end
grid on;
%     legend(legendstr);
xlabel('Time (s)');
ylabel('Absolute value');
title('Transmitted waveforms(PRS) of target(LTE BS)');

end