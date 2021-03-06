function [S] = gui_component_layout_sensor_fixed(S)

S.fh(2) = figure('units','pixels',...
    'position',[S.second_fig_pos_left 500 350 430],... % [1070 500 350 430]
    'name','input for sensor fixed',...
    'menubar','none',...
    'numbertitle','off',...
    'resize','off',...
    'closerequestfcn',{@fh_sub_crfcn, S.fh(1)});

S.tx(1) = uicontrol('units','pixels',...
    'style','text',...
    'unit','pix',...
    'position',[70 385 210 25],...
    'string','input for sensor fixed',...
    'fontweight','bold',...
    'fontsize',10, ...
    'foregroundcolor', 'b', ...
    'backgroundcolor',get(S.fh(2),'color'));

S.ed(1) = uicontrol('units','pixels',...
    'style','edit',...
    'unit','pix',...
    'string', '100', ...
    'position',[170 210 150 25]);
S.tx(2) = uicontrol('units','pixels',...
    'style','text',...
    'unit','pix',...
    'position',[30 205 110 25],...
    'string','trial length',...
    'fontweight','bold',...
    'backgroundcolor',get(S.fh(2),'color'));

S.ed(2) = uicontrol('units','pixels',...
    'style','edit',...
    'unit','pix',...
    'string', '1.5', ...
    'position',[170 245 150 25]);
S.tx(3) = uicontrol('units','pixels',...
    'style','text',...
    'unit','pix',...
    'position',[30 240 110 25],...
    'string','radius ratio',...
    'tooltipstring', 'distance ratio between target and farthest sensor from [0,0]', ...
    'fontweight','bold',...
    'backgroundcolor',get(S.fh(2),'color'));

S.ed(3) = uicontrol('units','pixels',...
    'style','edit',...
    'unit','pix',...
    'string', '100', ...
    'position',[170 280 150 25]);
S.tx(4) = uicontrol('units','pixels',...
    'style','text',...
    'unit','pix',...
    'position',[30 275 110 25],...
    'string','y of target position',...
    'tooltipstring', 'unit = meter', ...
    'fontweight','bold',...
    'backgroundcolor',get(S.fh(2),'color'));

S.ed(4) = uicontrol('units','pixels',...
    'style','edit',...
    'unit','pix',...
    'string', '100', ...
    'position',[170 315 150 25]);
S.tx(5) = uicontrol('units','pixels',...
    'style','text',...
    'unit','pix',...
    'position',[30 310 110 25],...
    'string','x of target position',...
    'tooltipstring', 'unit = meter', ...
    'fontweight','bold',...
    'backgroundcolor',get(S.fh(2),'color'));

S.ed(5) = uicontrol('units','pixels',...
    'style','edit',...
    'unit','pix',...
    'string', '5', ...
    'position',[170 350 150 25]);
S.tx(6) = uicontrol('units','pixels',...
    'style','text',...
    'unit','pix',...
    'position',[30 345 110 25],...
    'string','SNR in dB',...
    'fontweight','bold',...
    'backgroundcolor',get(S.fh(2),'color'));

S.pb(2) = uicontrol('style','push',...
    'unit','pix',...
    'position',[100 20 150 30],...
    'fontweight','bold',...
    'string','run sensor fixed',...
    'callback',{@pb_run_sensor_fixed_call, S.fh(1)});

end