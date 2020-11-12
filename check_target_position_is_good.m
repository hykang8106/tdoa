function [good_target] = check_target_position_is_good(target_position, target_radius_meter)

good_target = 1;

if target_position(1) > target_radius_meter || ...
        target_position(1) < -target_radius_meter || ...
        target_position(2) > target_radius_meter || ...
        target_position(2) < -target_radius_meter    
    good_target = 0;
end

% if target_position(1) > target_radius_meter || ...
%         target_position(1) < -target_radius_meter || ...
%         target_position(2) > target_radius_meter || ...
%         target_position(2) < -target_radius_meter
%     fprintf('###### target position is out of range\n');
%     return;
% end

end
