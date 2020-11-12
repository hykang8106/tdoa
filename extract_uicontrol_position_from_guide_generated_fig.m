function [] = extract_uicontrol_position_from_guide_generated_fig(fig_filename)
% extract uicontrol position from GUIDE generated fig file
% ### used to set uicontrol 'position' property when creating gui code programmatically
%
% ########################################################################
% #### uicontrol "unit" property MUST BE pixel
% #### if uicontrol "string" property is untouched, line 55 give error
% #### to cure this, set "string" property to what you want 
% ########################################################################
%
% [usage]
% extract_uicontrol_position_from_guide_generated_fig('tdoa_guide.fig')

% get figure handle
f = openfig(fig_filename, 'new', 'invisible');
% f = hgload(filename);
% f = findobj(f,'Type','figure');

% fp: figure property, structure
fp = get(f);
% figure property name
fp_name = fieldnames(fp);
% fp_value = struct2cell(fp);
fprintf('#### figure\n');
position = fp.Position

% get uicontrol handle(child handle)
uic = allchild(f);
% uicontrol length(how many uicontrol is included in figure?)
uic_length = length(uic);

% uicp: uicontrol property, structure array whose length is child_length
uicp = get(uic);
uic_length = length(uicp); % this must be same as length(uic)

% uicontrol property name
uicp_name = fieldnames(uicp(1)); % assumed that uicontrol have all same property name
% uicontrol property length
uicp_length = length(uicp_name);

uicp_all_value = cell(uicp_length, uic_length);
uicp_style_value = cell(uic_length, 1);
% 4 = [left, bottom, width, height]
uicp_position_value = zeros(uic_length, 4);

for n = 1 : uic_length
    % without output(x), "get" display uicontrol property(lengthy!)
    x = get(uic(n)); 
    
    style = uicp(n).Style;
    fprintf('\n##### style = %s\n', style);
    uicp_style_value{n} = style;
    
    if strcmp(style, 'text') || strcmp(style, 'checkbox') || strcmp(style, 'pushbutton')
        string = uicp(n).String;
        fprintf('string = %s\n', char(string));
%         fprintf('string = %s\n', string);
    end
    
    position = fix(uicp(n).Position)
    uicp_position_value(n, :) = position;
    
%     extent = fix(uicp(n).Extent) % nothing special
end

uicp_position_value;

uicp_name;
uicp_style_value;
uicp_all_value;

end

