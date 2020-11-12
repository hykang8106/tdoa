function [] = learn_msgbox(msg_text, msg_title, msg_icon, create_mode, font_size, text_color)
% change font size and color of text in msgbox
%
% https://www.mathworks.com/matlabcentral/newsreader/view_thread/73331
%
% [input]
% - msg_text:
% - msg_title:
% - msg_icon: 'none'(default) | 'error' | 'help' | 'warn'
% - create_mode: 'modal' | 'nonmodal'(default)
%
% [usage]
% learn_msgbox('test completed, end of winter, spring come', 'info', 'help', 'modal', 14, 'b')
% learn_msgbox('test completed, end of winter, spring come', 'info', [], 'modal', 14, 'b')

% % points pixels
% point_pixel = ...
%     [
%     8	7
%     9	8
%     10	9
%     11	10
%     12	11
%     13	12
%     14	13
%     ];

if ~isempty(msg_icon)
    h = msgbox({msg_text}, msg_title, msg_icon, create_mode);
else
    h = msgbox({msg_text}, msg_title, create_mode);
end

if isempty(msg_icon)
    icon_pixel_width = 0;
else
    % image object, CData = [ (50 by 50) uint8 array]
    icon_pixel_width = 50;
end

text_length = length(msg_text);

% ### 1: hardcored, not nice
% idx = find(point_pixel(:, 1) == font_size);
% figure_width = text_length * point_pixel(idx, 2) + icon_pixel_width; 
figure_width = text_length * (font_size - 1) + icon_pixel_width; 

% try to set msgbox width enough to display msg_text
set(h, 'units', 'pixels');
figure_position = get(h, 'position');
figure_position(3) = figure_width;
set(h, 'position', figure_position);

% get all graphic object
H = findall(h);

% obj_length = length(H);
% for n = 1 : obj_length
%     disp('##########');
%     get(H(n));
% end

% filter only text handle
text_h = unique(findobj(H, 'type', 'text'));

% set font_size and text_color of msg_text
set(text_h, 'FontSize', font_size, 'fontweight', 'bold', 'Color', text_color);

% try to set msg_text to center of msgbox
set(text_h, 'units', 'pixels');
text_position = get(text_h, 'position');
if isempty(msg_icon)
    text_position(1) = text_position(1) + 15; % ### 15: hardcored, not nice
    set(text_h, 'position', text_position);
end

% filter uicontrol(pushbutton) handle
pb_h = unique(findobj(H, 'type', 'uicontrol'));

% try to set pushbutton to center of msgbox
set(pb_h, 'units', 'pixels');
pb_position = get(pb_h, 'position');
pb_position(1) = figure_position(3) / 2.25; % ### 2.25: hardcored, not nice
set(pb_h, 'position', pb_position);

uiwait(h);

end

