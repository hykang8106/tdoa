function my_slider_using_appdata()

hfig = figure();
setappdata(hfig,'slidervalue',0);
setappdata(hfig,'difference',1);

slider = uicontrol('Parent', hfig,'Style','slider',...
    'Units','normalized',...
    'Position',[0.3 0.5 0.4 0.1],...
    'Tag','slider1',...
    'Callback',@slider_callback);

button = uicontrol('Parent', hfig,'Style','pushbutton',...
    'Units','normalized',...
    'Position',[0.4 0.3 0.2 0.1],...
    'String','Display Values',...
    'Callback',@button_callback);

end

function slider_callback(hObject,eventdata)

sval = get(hObject,'Value');
diffMax = get(hObject,'Max') - sval;
hfig = get(hObject,'Parent');
setappdata(hfig,'slidervalue',sval);
setappdata(hfig,'difference',diffMax);

end

function button_callback(hObject,eventdata)

hfig = get(hObject,'Parent');
currentval = getappdata(hfig,'slidervalue');
diffval = getappdata(hfig,'difference');
display([currentval diffval]);

end