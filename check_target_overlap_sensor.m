function [overlapped] = ...
    check_target_overlap_sensor(sensor_position, target_position, overlap_criterion_meter)

overlapped = 0;

sensor_length = size(sensor_position, 1);
for n = 1 : sensor_length
%     sensor_position(n, :);
%     sensor_position(n, :) - target_position;
    if sqrt(sum((sensor_position(n, :) - target_position).^2)) < overlap_criterion_meter
        overlapped = 1;
        break;
    end
end

end


