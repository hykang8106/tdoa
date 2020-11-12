function [S] = gui_component_layout_batch_bw_corr_snr(S)

S.fh(2) = figure('units','pixels',...
    'position',[S.second_fig_pos_left   290   361   510],... % [520   290   361   510]
    'name','input for batch bw corr snr',...
    'menubar','none',...
    'numbertitle','off',...
    'resize','off',...
    'closerequestfcn',{@fh_sub_crfcn, S.fh(1)});

S.tx(1) = uicontrol('units','pixels',...
    'style','text',...
    'unit','pix',...
    'position',[100   477   151    15],...
    'string','input for batch bw corr snr',...
    'fontweight','bold',...
    'fontsize',10, ...
    'foregroundcolor', 'b', ...
    'backgroundcolor',get(S.fh(2),'color'));

S.ed(1) = uicontrol('units','pixels',...
    'style','edit',...
    'unit','pix',...
    'string', '100, 100', ...
    'position',[200   437   101    23]);
S.tx(2) = uicontrol('units','pixels',...
    'style','text',...
    'unit','pix',...
    'position',[50   441   120    15],...
    'string','target position (x, y)',...
    'fontweight','bold',...
    'backgroundcolor',get(S.fh(2),'color'));

S.ed(2) = uicontrol('units','pixels',...
    'style','edit',...
    'unit','pix',...
    'string', '1.5', ...
    'position',[200   387   101    23]);
S.tx(3) = uicontrol('units','pixels',...
    'style','text',...
    'unit','pix',...
    'position',[50   391   101    15],...
    'string','radius ratio',...
    'tooltipstring', 'distance ratio between target and farthest sensor from [0,0]', ...
    'fontweight','bold',...
    'backgroundcolor',get(S.fh(2),'color'));

S.ed(3) = uicontrol('units','pixels',...
    'style','edit',...
    'unit','pix',...
    'string', '100', ...
    'position',[200   337   101    23]);
S.tx(4) = uicontrol('units','pixels',...
    'style','text',...
    'unit','pix',...
    'position',[50   341   101    15],...
    'string','trial length',...
    'fontweight','bold',...
    'backgroundcolor',get(S.fh(2),'color'));

S.lb(1) = uicontrol('units','pixels',...
    'style','listbox',...
    'unit','pix',...
    'string', {'6', '15', '25'}, ...
    'value', 2, ...
    'position',[200   242   101    68]);
S.tx(5) = uicontrol('units','pixels',...
    'style','text',...
    'unit','pix',...
    'position',[50   273   101    15],...
    'string','NDLRB',...
    'fontweight','bold',...
    'backgroundcolor',get(S.fh(2),'color'));

S.ed(4) = uicontrol('units','pixels',...
    'style','edit',...
    'unit','pix',...
    'string', '-5, 0, 5', ...
    'position',[200   187   101    23]);
S.tx(6) = uicontrol('units','pixels',...
    'style','text',...
    'unit','pix',...
    'position',[50   191   101    15],...
    'string','SNR in dB (vector)',...
    'fontweight','bold',...
    'backgroundcolor',get(S.fh(2),'color'));

S.ed(5) = uicontrol('units','pixels',...
    'style','edit',...
    'unit','pix',...
    'string', '1:15', ...
    'position',[200   137   101    23]);
S.tx(7) = uicontrol('units','pixels',...
    'style','text',...
    'unit','pix',...
    'position',[50   141   101    15],...
    'string','NPRSRB (vector)',...
    'fontweight','bold',...
    'backgroundcolor',get(S.fh(2),'color'));

S.ed(6) = uicontrol('units','pixels',...
    'style','edit',...
    'unit','pix',...
    'string', '1, 2', ...
    'position',[200    87   101    23]);
S.tx(8) = uicontrol('units','pixels',...
    'style','text',...
    'unit','pix',...
    'position',[50    91   101    15],...
    'string','nsubframe (vector)',...
    'fontweight','bold',...
    'backgroundcolor',get(S.fh(2),'color'));

S.pb(2) = uicontrol('style','push',...
    'unit','pix',...
    'position',[100    37   151    23],...
    'fontweight','bold',...
    'string','run batch bw corr snr',...
    'callback',{@pb_run_batch_bw_corr_snr_call, S.fh(1)});

end
