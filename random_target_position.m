function [target_position] = random_target_position(uca_radius_meter, radius_ratio)

r = rand * uca_radius_meter * radius_ratio;
t = rand * 2 * pi;
target_position = [r * cos(t), r * sin(t)];

% 
% True_position = zeros(Trials, 2);
% Est_position = zeros(Trials, 2);
% 
% % Generate position of source
% for i = 1 : Trials
%     r = rand * Radius * radius_ratio;
%     t = rand * 2 * pi;
%     True_position(i, 1) = r * cos(t);
%     True_position(i, 2) = r * sin(t);
% end
% True_position;

end