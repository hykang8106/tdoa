function my_slider_using_guidata()
% learn guidata

hfig = figure();
% data = struct('val',0,'diffMax',1,'xyz',[]);
% guidata(hfig, data);
slider = uicontrol('Parent', hfig,'Style','slider',...
    'Units','normalized',...
    'Position',[0.3 0.5 0.4 0.1],...
    'Tag','slider1',...
    'Callback',@sub_slider_callback_using_guidata);

button = uicontrol('Parent', hfig,'Style','pushbutton',...
    'Units','normalized',...
    'Position',[0.4 0.3 0.2 0.1],...
    'String','Display Values',...
    'Tag','button1',...
    'Callback',@sub_button_callback_using_guidata);

data = struct('val',0,'diffMax',1,'xyz',[]);
guidata(hfig, data);

guihandles(button)

end

% function sub_slider_callback_using_guidata(varargin)
% % varargin = {hObject, eventdata}
% 
% hObject = varargin{1};
% data = guidata(hObject);
% data.val = get(hObject,'Value');
% data.diffMax = get(hObject,'Max') - data.val;
% data.xyz = 10;
% guidata(hObject,data);
% 
% end

% function sub_button_callback_using_guidata(varargin)
% 
% hObject = varargin{1};
% data = guidata(hObject);
% display([data.val data.diffMax]);
% display(data.xyz);
% 
% end

