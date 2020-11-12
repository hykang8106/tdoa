function sub_slider_callback_using_guidata(varargin)
% #### used in my_slider_using_guidata.m
%
% varargin = {hObject, eventdata}

hObject = varargin{1};
data = guidata(hObject);
data.val = get(hObject,'Value');
data.diffMax = get(hObject,'Max') - data.val;
data.xyz = 10;
guidata(hObject,data);

end

