function [ref_sensor] = choose_reference_sensor(rx)
% choose reference sensor which is nearest to target emitter, so largest signal is received 

if iscell(rx)
    rx = cell2mat(rx);
end
% whos rx_array

% figure;
% plot(abs(rx_array));

[Y, ref_sensor] = max(sum(abs(rx)));
% sum(abs(rx));

end

