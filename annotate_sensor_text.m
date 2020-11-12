function [] = annotate_sensor_text(x, y, sensor_position)
% #### called in plot_tdoa_result_sensor_fixed.m to annotate sensor text on plot

sensor_length = size(sensor_position, 1);

x_length = length(x);
y_length = length(y);

% [X, Y] = meshgrid(x, y);

for n = 1 : sensor_length
    sp = sensor_position(n, :);

    % #### comment out meshgrid: 
    % #### use simple method, just put sensor name text on sensor position
    
%     % ## repmat input(y_length, x_length) is right order? (y_length, x_length) vs (x_length, y_length)
%     Z = (X - repmat(sp(1), y_length, x_length)).^2 + (Y - repmat(sp(2), y_length, x_length)).^2;
%     min_Z = min(min(Z));
%     
%     [I, J] = find(Z == min_Z);
%     
%     x_position = x(J);
%     y_position = y(I);

    x_position = sp(1);
    y_position = sp(2);
    
    h = text(x_position, y_position, sprintf('s%d', n), ...
        'VerticalAlignment', 'bottom', 'color', 'r', 'fontweight', 'bold');
    
    % ##### bring 'sensor' text to front
    % ##### reference: https://kr.mathworks.com/matlabcentral/newsreader/view_thread/292321
    set(h,'erasemode','xor');
    set(h,'erasemode','background');
    
end

end