function [uicp_position_value] = extract_uicontrol_position_from_guide_generated_fig(fig_filename)
% extract uicontrol position from GUIDE generated fig file
% ### when creating gui code programmatically, used to set uicontrol 'position' property
%
% [usage]
% extract_uicontrol_position_from_guide_generated_fig('tdoa_guide.fig')

% get figure handle
f = openfig(fig_filename, 'reuse', 'invisible');
% f = hgload(filename);
% f = findobj(f,'Type','figure');

% fp: figure property, structure
fp = get(f);
% figure property name
fp_name = fieldnames(fp);
% fp_value = struct2cell(fp);
fp.Position

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
uicp_position_value = zeros(uic_length, 4);

for n = 1 : uic_length
    X = get(uic(n));
    uicp_style_value{n} = uicp(n).Style;
    
    uicp_position_value(n, :) = uicp(n).Position;
    
    % set empty to callback function property
    uicp(n).ButtonDownFcn = '';
    uicp(n).Callback = '';
    uicp(n).CreateFcn = '';
    uicp(n).DeleteFcn = '';
    uicp(n).KeyPressFcn = '';
    
    uicp(n).BackgroundColor = num2str(uicp(n).BackgroundColor);
    uicp(n).Extent = num2str(uicp(n).Extent);
    uicp(n).ForegroundColor = num2str(uicp(n).ForegroundColor);
    uicp(n).Position = num2str(uicp(n).Position);
    uicp(n).SliderStep = num2str(uicp(n).SliderStep);
    
    % ###### caution: NOT {}, BUT () #######
    uicp_all_value(:, n) = struct2cell(uicp(n));
end

uicp_position_value

uicp_name;
uicp_style_value;
uicp_all_value;

% [p, name, e] = fileparts(fig_filename);
% excel_filename = sprintf('%s.xlsx', name);
% 
% xlswrite(excel_filename, uicontrol_property_name, 1, 'a5');
% xlswrite(excel_filename, uicontrol_property_value, 1, 'b5');
% xlswrite(excel_filename, uicontrol_property_style_value, 1, 'b4');

end

