function sub_button_callback_using_guidata(varargin)

hObject = varargin{1};
data = guidata(hObject);
display([data.val data.diffMax]);
display(data.xyz);

end