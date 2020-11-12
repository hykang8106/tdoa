function my_slider_using_nested_callback()

hfig = figure();
data = struct('val',0,'diffMax',1);
slider = uicontrol('Parent', hfig,'Style','slider',...
    'Units','normalized',...
    'Position',[0.3 0.5 0.4 0.1],...
    'Tag','slider1',...
    'Callback',@slider_callback);

button = uicontrol('Parent',hfig,'Style','pushbutton',...
    'Units','normalized',...
    'Position',[0.4 0.3 0.2 0.1],...
    'String','Display Difference',...
    'Callback',@button_callback);

    function slider_callback(hObject,eventdata)
        
        sval = get(hObject,'Value');
        diffMax = get(hObject,'Max') - sval;
        data.val = sval;
        data.diffMax = diffMax;
        
    end

    function button_callback(hObject,eventdata)
        
        display([data.val data.diffMax]);
        
    end
end