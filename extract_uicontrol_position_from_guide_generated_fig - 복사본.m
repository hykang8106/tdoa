function [] = extract_uicontrol_position_from_guide_generated_fig(fig_filename)
% extract uicontrol position from GUIDE generated fig file
% ### when creating gui code programmatically, used to set uicontrol 'position' property
%
% [usage]
% extract_uicontrol_position_from_guide_generated_fig('tdoa_guide.fig')

f = openfig(fig_filename, 'reuse', 'invisible');
% f = hgload(filename);
% f = findobj(f,'Type','figure');

% FP: structure
FP = get(f);
figure_property_name = fieldnames(FP);
figure_property_value = struct2cell(FP);
FP.Position

c = allchild(f);
child_length = length(c);

% CP: structure array whose length is child_length
CP = get(c);
child_length = length(CP); % this must be same as length(c)

uicontrol_property_name = fieldnames(CP);
uicontrol_property_length = length(uicontrol_property_name);
uicontrol_property_value = cell(uicontrol_property_length, child_length);
uicontrol_property_style_value = cell(1, child_length);

uicontrol_position = cell(child_length, 5);

for n = 1 : child_length
    uicontrol_property_style_value{n} = CP(n).Style;
    X = get(c(n));
    
    uicontrol_position{n, 1} = CP(n).Style;
    uicontrol_position{n, 2} = CP(n).Position(1);
    uicontrol_position{n, 3} = CP(n).Position(2);
    uicontrol_position{n, 4} = CP(n).Position(3);
    uicontrol_position{n, 5} = CP(n).Position(4);
    
    % set empty to callback function property
    CP(n).ButtonDownFcn = '';
    CP(n).Callback = '';
    CP(n).CreateFcn = '';
    CP(n).DeleteFcn = '';
    CP(n).KeyPressFcn = '';
    
    CP(n).BackgroundColor = num2str(CP(n).BackgroundColor);
    CP(n).Extent = num2str(CP(n).Extent);
    CP(n).ForegroundColor = num2str(CP(n).ForegroundColor);
    CP(n).Position = num2str(CP(n).Position);
    CP(n).SliderStep = num2str(CP(n).SliderStep);
    
    % ###### caution: NOT {}, BUT () #######
    uicontrol_property_value(:, n) = struct2cell(CP(n));
end

uicontrol_position

uicontrol_property_name;
uicontrol_property_style_value;
uicontrol_property_value;

% [p, name, e] = fileparts(fig_filename);
% excel_filename = sprintf('%s.xlsx', name);
% 
% xlswrite(excel_filename, uicontrol_property_name, 1, 'a5');
% xlswrite(excel_filename, uicontrol_property_value, 1, 'b5');
% xlswrite(excel_filename, uicontrol_property_style_value, 1, 'b4');

end

