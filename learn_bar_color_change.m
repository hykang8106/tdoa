function [] = learn_bar_color_change(steps, pause_sec)
% learn wait bar color change
%
% [usage]
% learn_bar_color_change(1000, .01)
%

if nargin ~= 2
    fprintf(2, '#### [usage] learn_bar_color_change(steps, pause_sec)\n');
    return;
end

rgb_map = [
    0 0 1 % 'b'
    0 1 0 % 'g'
    1 0 0 % 'r'
    0 1 1 % 'c'
    1 0 1 % 'm'
    1 1 0 % 'y'
    0 0 0 % 'k'
    1 1 1 % 'w'
    ];

rgb_map = jet(512);
% rgb_vec = bone(512);
rgb_len = size(rgb_map, 1);

h = waitbar(0, 'Please wait...', 'name', 'changing bar color');
hw = findobj(h, 'Type', 'Patch');

% edge_color = [0 1 0];
% face_color = [0 1 0];
% 
% set(hw, 'EdgeColor', edge_color, 'FaceColor', face_color); % changes the color to green
% % set(hw, 'EdgeColor', [0 1 0], 'FaceColor', [0 1 0]); % changes the color to green

u = tic;
for i = 1 : steps
    w = toc(u);
    % reference:
    % http://stackoverflow.com/questions/12210583/is-there-a-matlab-function-to-convert-elapsed-seconds-to-hhmmss-format
    z = fix(mod(w, [0, 3600, 60])./[3600, 60, 1]);
    
    waitbar(i / steps, h, sprintf('%.0f %%, elapsed time = %d : %02d : %02d', ...
        i / steps * 100, z(1), z(2), z(3)));
    
    rgb_idx = mod(i, rgb_len);
    if ~rgb_idx
        rgb_idx = rgb_len;
    end
    edge_color = rgb_map(rgb_idx, :);
    face_color = rgb_map(rgb_idx, :);
    
    % reference:
    % https://kr.mathworks.com/matlabcentral/answers/101667-how-do-i-change-the-color-of-the-bar-in-a-waitbar-in-matlab-7-5-r2007b
    set(hw, 'EdgeColor', edge_color, 'FaceColor', face_color);
    
    pause(pause_sec);
end

close(h);

end



