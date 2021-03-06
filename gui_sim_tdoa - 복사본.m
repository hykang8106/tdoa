function [] = gui_sim_tdoa()
% ####################################################
% ##### gui version of simulate tdoa location
% ##### reference code: GUI_41.m (Author:  Matt Fig)
% ####################################################
%

% radio_button_length = 8;
% figure_length = radio_button_length + 1;

S.sensor_position = [];

S.rician_param = [];

S.uca_radius_meter = [];

S.fading_param_filename = [];

S.str_ndlrb = {'6';'15';'25'};  % String for popups.

% #############################################################################
% ### temporary: must rewrite,
% ### nprsrb is dependant on ndlrb, nprsrb = 1 : current value of ndlrb
% #############################################################################
S.str_nprsrb = {'1';'2';'3';'4';'5';'6'};  % String for popups.

S.str_subframe_length = {'1';'2'};  % String for popups.

% S.X = -10:.01:10;  % The X values for plotting.
S.fh(1) = figure('units','pixels',...
    'position',[200 650 850 250],...
    'menubar','none',...
    'name','simulation of tdoa location',...
    'numbertitle','off',...
    'resize','off',...
    'closerequestfcn',{@fh_main_crfcn});
S.txmain = uicontrol('units','pixels',...
    'style','text',...
    'unit','pix',...
    'position',[200 190 450 25],...
    'string','Radio Emitter Location using TDOA Sensors',...
    'fontweight','bold',...
    'backgroundcolor',get(S.fh(1),'color'), ...
    'fontsize',15, ...
    'foregroundcolor', 'b');
S.bg = uibuttongroup('units','pix',...
    'pos',[20 70 810 90]);
S.rd(1) = uicontrol(S.bg,...
    'style','rad',...
    'unit','pix',...
    'position',[20 50 130 30],...
    'fontweight','bold',...
    'string','sensor fixed');
S.SEL = 1;  % The selectedobject property of S.bg
S.SEL_old = 1;
S.rd(2) = uicontrol(S.bg,...
    'style','rad',...
    'unit','pix',...
    'position',[20 10 130 30],...
    'fontweight','bold',...
    'string','target fixed');
S.rd(3) = uicontrol(S.bg,...
    'style','rad',...
    'unit','pix',...
    'position',[220 50 130 30],...
    'fontweight','bold',...
    'string','batch sensor fixed');
S.rd(4) = uicontrol(S.bg,...
    'style','rad',...
    'unit','pix',...
    'position',[220 10 130 30],...
    'fontweight','bold',...
    'string','batch target fixed');
S.rd(5) = uicontrol(S.bg,...
    'style','rad',...
    'unit','pix',...
    'position',[420 50 130 30],...
    'fontweight','bold',...
    'string','batch bw corr snr');
S.rd(6) = uicontrol(S.bg,...
    'style','rad',...
    'unit','pix',...
    'position',[420 10 130 30],...
    'fontweight','bold',...
    'string','plot target fixed');
S.rd(7) = uicontrol(S.bg,...
    'style','rad',...
    'unit','pix',...
    'position',[620 50 130 30],...
    'fontweight','bold',...
    'string','plot sensor fixed');
S.rd(8) = uicontrol(S.bg,...
    'style','rad',...
    'unit','pix',...
    'position',[620 10 130 30],...
    'fontweight','bold',...
    'string','plot bw corr snr');
S.pb_main = uicontrol('style','push',...
    'unit','pix',...
    'position',[355 20 150 30],...
    'fontweight','bold',...
    'string','input parameters',...
    'callback',{@pb_input_param_call});
S.fm(1) = uimenu(S.fh,...
    'label','load sensor position',...
    'callback',{@fm_call},...
    'enable','on');
S.fm(2) = uimenu(S.fh,...
    'label','load fading parameter',...
    'callback',{@fm_call},...
    'enable','on');

    function [] = fm_call(varargin)
        % Callback for the figure menu.
        switch gcbo
            case S.fm(1)
                S.sensor_position_filename = uigetfile('*.txt','select sensor position file');
                try
                    [S.sensor_position, S.uca_radius_meter] = ...
                        get_sensor_position_from_file(S.sensor_position_filename);
                catch
                    disp('Unable to Load.  Check Name and Try Again.')
                end
            case S.fm(2)
                S.fading_param_filename = uigetfile('*.xlsx','select fading parameter file');
                try
                    % ############### not nice, rewrite ##########
                    xlsrange = 'b3:d7';
                    S.rician_param = ...
                        get_rician_parameter_from_excel_file(S.fading_param_filename, xlsrange);
                    if size(S.rician_param, 1) ~= size(S.sensor_position, 1)
                        error('##### row length of rician parameter must be same as sensor length\n');
                    end
                catch
                    disp('Unable to Load.  Check Name and Try Again.')
                end
            otherwise
        end
    end

    function [] = pb_input_param_call(varargin)
        % Callback for the pushbutton.
        sel = findobj(get(S.bg,'selectedobject'));  % See BUG note in GUI_8
        sel;
        S.fh;
        S.SEL_old;
        S.SEL = find(S.rd == sel);  % Store current radiobutton.
        if (length(S.fh) == 2) && (S.SEL ~= S.SEL_old)
            delete(S.fh(2));
            S.fh(2) = [];
        end
        S.SEL_old = S.SEL;
        S.SEL;
        
        switch sel
            
            % ##########################
            % ##### sensor fixed
            % ##########################
            case S.rd(1) 
                %                 fprintf('fig length = %d\n', length(S.fh));
                if length(S.fh) == 1  % We haven't been here before.

                    % ####### failed to make this code function
                    % ####### error message: can't find callback function
                    % ####### why?
%                     gui_populate_sensor_fixed_uicontrol(S);

                    % This is where we make the other two figures.
                    S.fh(2) = figure('units','pixels',...
                        'position',[1070 500 350 430],...
                        'name','input for sensor fixed',...
                        'menubar','none',...
                        'numbertitle','off',...
                        'resize','off',...
                        'closerequestfcn',{@fh_sub_crfcn});
                    
                    S.tx = uicontrol('units','pixels',...
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
                        'tooltipstring', 'ratio between target distance and farthest sensor distance from [0,0]', ...
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
                    
                    S.pp(1) = uicontrol('style','pop',...
                        'unit','pix',...
                        'position',[170 80 70 30],...
                        'string',S.str_subframe_length, ...
                        'value',1);
                    S.tx(7) = uicontrol('units','pixels',...
                        'style','text',...
                        'unit','pix',...
                        'position',[30 80 110 25],...
                        'string','subframe length',...
                        'tooltipstring', 'determine correlation length', ...
                        'fontweight','bold',...
                        'backgroundcolor',get(S.fh(2),'color'));
                    
                    %                     ca = get(S.pp(3), {'string', 'value'});
                    %                     ndlrb = str2num(ca{1}{ca{2}})
                    %                     cell_nprsrb = cell(ndlrb, 1);
                    %                     for n = 1 : ndlrb
                    %                         cell_nprsrb{n} = num2str(n);
                    %                     end
                    %                     S.str_nprsrb = cell_nprsrb;
                    S.pp(2) = uicontrol('style','pop',...
                        'unit','pix',...
                        'position',[170 120 70 30],...
                        'string',S.str_nprsrb,...
                        'value',1);
                    S.tx(8) = uicontrol('units','pixels',...
                        'style','text',...
                        'unit','pix',...
                        'position',[30 120 110 25],...
                        'string','NPRSRB',...
                        'tooltipstring', 'Number of PRS Resource Block, determine signal bandwidth', ...
                        'foregroundcolor', 'r', ...
                        'fontweight','bold',...
                        'backgroundcolor',get(S.fh(2),'color'));
                    
                    S.pp(3) = uicontrol('style','pop',...
                        'unit','pix',...
                        'position',[170 160 70 30],...
                        'string',S.str_ndlrb,...
                        'value',1);
                    S.tx(9) = uicontrol('units','pixels',...
                        'style','text',...
                        'unit','pix',...
                        'position',[30 160 110 25],...
                        'string','NDLRB',...
                        'tooltipstring', 'Number of DownLink Resource Block, determine signal bandwidth', ...
                        'fontweight','bold',...
                        'backgroundcolor',get(S.fh(2),'color'));
                    
                    S.pb_sensor_fixed = uicontrol('style','push',...
                        'unit','pix',...
                        'position',[100 20 150 30],...
                        'fontweight','bold',...
                        'string','run sensor fixed',...
                        'callback',{@pb_sensor_fixed_call});
                    
                end
                
            case S.rd(2)  % target fixed
                %                 if length(S.fh) == 1  % We haven't been here before.
                %                     % This is where we make the other two figures.
                %                     S.fh(2) = figure('units','pixels',...
                %                         'position',[330 170 350 430],...
                %                         'name','input for sensor fixed',...
                %                         'menubar','none',...
                %                         'numbertitle','off',...
                %                         'resize','off',...
                %                         'closerequestfcn',{@fh_sub_crfcn});
                %                 end
                %                 disp('not implemented');
            
            % ##########################
            % ##### batch sensor fixed
            % ##########################
            case S.rd(3)
                if length(S.fh) == 1  % We haven't been here before.
                    % This is where we make the other two figures.
                    S.fh(2) = figure('units','pixels',...
                        'position',[1070 500 350 430],...
                        'name','input for batch sensor fixed',...
                        'menubar','none',...
                        'numbertitle','off',...
                        'resize','off',...
                        'closerequestfcn',{@fh_sub_crfcn});
                    
                    S.tx = uicontrol('units','pixels',...
                        'style','text',...
                        'unit','pix',...
                        'position',[70 385 210 25],...
                        'string','input for batch sensor fixed',...
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
                        'tooltipstring', 'ratio between target distance and farthest sensor distance from [0,0]', ...
                        'fontweight','bold',...
                        'backgroundcolor',get(S.fh(2),'color'));
                    
%                     S.ed(3) = uicontrol('units','pixels',...
%                         'style','edit',...
%                         'unit','pix',...
%                         'string', '100', ...
%                         'position',[170 280 150 25]);, 
%                     S.tx(4) = uicontrol('units','pixels',...
%                         'style','text',...
%                         'unit','pix',...
%                         'position',[30 275 110 25],...
%                         'string','y of target position',...
%                         'fontweight','bold',...
%                         'backgroundcolor',get(S.fh(2),'color'));
                    
                    S.ed(4) = uicontrol('units','pixels',...
                        'style','edit',...
                        'unit','pix',...
                        'string', '100', ...
                        'position',[170 315 150 25]);
                    S.tx(5) = uicontrol('units','pixels',...
                        'style','text',...
                        'unit','pix',...
                        'position',[30 310 110 25],...
                        'string','target position step',...
                        'tooltipstring', 'same step is applied to x and y, unit = meter', ...
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
                    
                    S.pp(1) = uicontrol('style','pop',...
                        'unit','pix',...
                        'position',[170 80 70 30],...
                        'string',S.str_subframe_length, ...
                        'value',1);
                    S.tx(7) = uicontrol('units','pixels',...
                        'style','text',...
                        'unit','pix',...
                        'position',[30 80 110 25],...
                        'string','subframe length',...
                        'fontweight','bold',...
                        'backgroundcolor',get(S.fh(2),'color'));
                    
                    %                     ca = get(S.pp(3), {'string', 'value'});
                    %                     ndlrb = str2num(ca{1}{ca{2}})
                    %                     cell_nprsrb = cell(ndlrb, 1);
                    %                     for n = 1 : ndlrb
                    %                         cell_nprsrb{n} = num2str(n);
                    %                     end
                    %                     S.str_nprsrb = cell_nprsrb;
                    S.pp(2) = uicontrol('style','pop',...
                        'unit','pix',...
                        'position',[170 120 70 30],...
                        'string',S.str_nprsrb,...
                        'value',1);
                    S.tx(8) = uicontrol('units','pixels',...
                        'style','text',...
                        'unit','pix',...
                        'position',[30 120 110 25],...
                        'string','NPRSRB',...
                        'tooltipstring', 'Number of PRS Resource Block', ...
                        'foregroundcolor', 'r', ...
                        'fontweight','bold',...
                        'backgroundcolor',get(S.fh(2),'color'));
                    
                    S.pp(3) = uicontrol('style','pop',...
                        'unit','pix',...
                        'position',[170 160 70 30],...
                        'string',S.str_ndlrb,...
                        'value',1);
                    S.tx(9) = uicontrol('units','pixels',...
                        'style','text',...
                        'unit','pix',...
                        'position',[30 160 110 25],...
                        'string','NDLRB',...
                        'tooltipstring', 'Number of DownLink Resource Block', ...
                        'fontweight','bold',...
                        'backgroundcolor',get(S.fh(2),'color'));
                    
                    S.pb_batch_sensor_fixed = uicontrol('style','push',...
                        'unit','pix',...
                        'position',[100 20 150 30],...
                        'fontweight','bold',...
                        'string','run batch sensor fixed',...
                        'callback',{@pb_batch_sensor_fixed_call});
                    
                end
                
            case S.rd(4)  % batch target fixed
                %                 if length(S.fh) == 1  % We haven't been here before.
                %                     % This is where we make the other two figures.
                %                     S.fh(2) = figure('units','pixels',...
                %                         'position',[330 170 350 430],...
                %                         'name','input for sensor fixed',...
                %                         'menubar','none',...
                %                         'numbertitle','off',...
                %                         'resize','off',...
                %                         'closerequestfcn',{@fh_sub_crfcn});
                %                 end
                %                 S.COL = 'm';
            case S.rd(5)  % batch bw corr snr
                %                 if length(S.fh) == 1  % We haven't been here before.
                %                     % This is where we make the other two figures.
                %                     S.fh(2) = figure('units','pixels',...
                %                         'position',[330 170 350 430],...
                %                         'name','input for sensor fixed',...
                %                         'menubar','none',...
                %                         'numbertitle','off',...
                %                         'resize','off',...
                %                         'closerequestfcn',{@fh_sub_crfcn});
                %                 end
                %                 S.COL = 'm';
            case S.rd(6)  % plot target fixed
                %                 if length(S.fh) == 1  % We haven't been here before.
                %                     % This is where we make the other two figures.
                %                     S.fh(2) = figure('units','pixels',...
                %                         'position',[330 170 350 430],...
                %                         'name','input for sensor fixed',...
                %                         'menubar','none',...
                %                         'numbertitle','off',...
                %                         'resize','off',...
                %                         'closerequestfcn',{@fh_sub_crfcn});
                %                 end
                %                 S.COL = 'm';
            case S.rd(7)  % plot sensor fixed
                %                 if length(S.fh) == 1  % We haven't been here before.
                %                     % This is where we make the other two figures.
                %                     S.fh(2) = figure('units','pixels',...
                %                         'position',[330 170 350 430],...
                %                         'name','input for sensor fixed',...
                %                         'menubar','none',...
                %                         'numbertitle','off',...
                %                         'resize','off',...
                %                         'closerequestfcn',{@fh_sub_crfcn});
                %                 end
                %                 S.COL = 'm';
            case S.rd(8)  % plot bw corr snr
                %                 if length(S.fh) == 1  % We haven't been here before.
                %                     % This is where we make the other two figures.
                %                     S.fh(2) = figure('units','pixels',...
                %                         'position',[330 170 350 430],...
                %                         'name','input for sensor fixed',...
                %                         'menubar','none',...
                %                         'numbertitle','off',...
                %                         'resize','off',...
                %                         'closerequestfcn',{@fh_sub_crfcn});
                %                 end
                %                 S.COL = 'm';
            otherwise
                % Very unlikely I think.
        end
        
        function [] = pb_batch_sensor_fixed_call(varargin)
            
            % ### tried to make this code function(gui_get_sensor_position_dialog_box.m),
            % ### but when go to "catch" section, sensor_position is still empty, this cause error
            if isempty(S.sensor_position)
                uiwait(msgbox({'empty sensor position', 'dialog box will be opened'},...
                    'warning','modal'));
                sensor_position_filename = uigetfile('*.txt','select sensor position file');
                try
                    [S.sensor_position, S.uca_radius_meter] = ...
                        get_sensor_position_from_file(sensor_position_filename);
                catch
                    disp('Unable to Load.  Check Name and Try Again.');
                    return;
                end
            end
            
            % disable load sensor position file and fading parameter file
            set(S.fm,{'enable'},{'off';'off'});
            
            trial_length = str2num(get(S.ed(1),'String')); % trial length
            % get 'string', 'value' property of edit uicontrol
            get(S.ed(1), {'string', 'value'});
            % get output example = '100' [0], ##### what is [0]? [ans] no meaning in edit uicontrol, index begin from 1
            
            radius_ratio = str2num(get(S.ed(2),'String')); % radius ratio
            
%             target_pos_y = str2num(get(S.ed(3),'String')) % y of target position
            
            target_pos_step = str2num(get(S.ed(4),'String')) % x of target position
            
            snr_db = str2num(get(S.ed(5),'String')); % snr in db
            
            % ca = cell array
            ca = get(S.pp(1), {'string', 'value'}); % subframe length
            subframe_length = str2num(ca{1}{ca{2}});
            
            ca = get(S.pp(2), {'string', 'value'}); % nprsrb
            nprsrb = str2num(ca{1}{ca{2}});
            
            ca = get(S.pp(3), {'string', 'value'}); % ndlrb
            ndlrb = str2num(ca{1}{ca{2}});
            
            
            disp('not implemented');
            
            
            % enable load sensor position file and fading parameter file
            set(S.fm,{'enable'},{'on';'on'});
            
        end
        
        function [] = pb_sensor_fixed_call(varargin)
            
            % ### tried to make this code function(gui_get_sensor_position_dialog_box.m),
            % ### but when go to "catch" section, sensor_position is still empty
            % ### this cause error
            if isempty(S.sensor_position)
                uiwait(msgbox({'empty sensor position', 'dialog box will be opened'},...
                    'warning','modal'));
                sensor_position_filename = uigetfile('*.txt','select sensor position file');
                try
                    [S.sensor_position, S.uca_radius_meter] = ...
                        get_sensor_position_from_file(sensor_position_filename);
                catch
                    disp('Unable to Load.  Check Name and Try Again.');
                    return;
                end
            end
            
            % disable load sensor position file and fading parameter file
            set(S.fm,{'enable'},{'off';'off'});
            
            [trial_length, radius_ratio, target_pos_y, target_pos_x, ...
            snr_db, subframe_length, nprsrb, ndlrb] = gui_sensor_fixed_input_from_uicontrol(S);
            
%             trial_length = str2num(get(S.ed(1),'String')); % trial length
%             % get 'string', 'value' property of edit uicontrol
%             get(S.ed(1), {'string', 'value'});
%             % get output example = '100' [0], ##### what is [0]? [ans] no meaning in edit uicontrol, index begin from 1
%             
%             radius_ratio = str2num(get(S.ed(2),'String')); % radius ratio
%             
%             target_pos_y = str2num(get(S.ed(3),'String')); % y of target position
%             
%             target_pos_x = str2num(get(S.ed(4),'String')); % x of target position
%             
%             snr_db = str2num(get(S.ed(5),'String')); % snr in db
%             
%             % ca = cell array
%             ca = get(S.pp(1), {'string', 'value'}); % subframe length
%             subframe_length = str2num(ca{1}{ca{2}});
%             
%             ca = get(S.pp(2), {'string', 'value'}); % nprsrb
%             nprsrb = str2num(ca{1}{ca{2}});
%             
%             ca = get(S.pp(3), {'string', 'value'}); % ndlrb
%             ndlrb = str2num(ca{1}{ca{2}});
            
            [bw_mhz, fs, nfft, sample_length] = get_bw_from_prs_spec_db(ndlrb, nprsrb, subframe_length);
            
            [tx_signal, fs, nfft] = generate_target_signal_lte_prs(ndlrb, nprsrb, subframe_length);
            
            plot_signal = 0; plot_position = 0; use_only_torrieri_method = 1;
            
            target_position = [target_pos_x, target_pos_y];
            
            target_radius_meter = S.uca_radius_meter * radius_ratio;
            
            position_error_torrieri = zeros(1, trial_length);
            cep = zeros(1, trial_length);
            gdop = zeros(1, trial_length);
            lambda1 = zeros(1, trial_length);
            lambda2 = zeros(1, trial_length);
            theta_rad = zeros(1, trial_length);
            % ###### for computing position estimation covariance
            x_est_torrieri = zeros(2, trial_length);
            
            for n = 1 : trial_length
                [position_error_torrieri(n), x_est_torrieri(:, n), cep(n), gdop(n), ...
                    lambda1(n), lambda2(n), theta_rad(n)] = ...
                    sub_simulate_tdoa_sensor_fixed(S.sensor_position, snr_db, target_position, ...
                    plot_signal, plot_position, use_only_torrieri_method, target_radius_meter, S.rician_param, ...
                    tx_signal, fs, nfft, bw_mhz);
            end
            
            % ###### target location beyond target_position_limit is ignored:
            % ###### excluded in histogram
            % ###### excluded in computing mean position error
            
            % #############################################################################################
            % ###### current target_position_limit may be too small
            % ###### target_position_limit = uca_radius_meter * 3 is good?
            % ###### target_position_limit also exist in batch_simulate_tdoa_bw_corr_snr.m
            % #############################################################################################
            target_position_limit = target_radius_meter;
            idx = position_error_torrieri <= target_position_limit;
            position_error_excl = position_error_torrieri(idx);
            
            mean_position_error = mean(position_error_excl);
            
            if isempty(S.rician_param)
                title_text = sprintf('[location] snr = %d db, trial = %d, location error = %d m, excl = %d', ...
                    snr_db, trial_length, round(mean_position_error), trial_length - sum(idx));
            else
                title_text = sprintf('[location, rician] snr = %d db, trial = %d, location error = %d m, excl = %d', ...
                    snr_db, trial_length, round(mean_position_error), trial_length - sum(idx));
            end
            plot_sensor_position_only(S.sensor_position, target_radius_meter, title_text, ...
                S.fading_param_filename);
            hold on;
            plot_target_position_and_ellipse(target_position, x_est_torrieri, ...
                cep, gdop, lambda1, lambda2, theta_rad);
            
            histogram_bin_length = 100;
            plot_position_error_histogram(position_error_excl, histogram_bin_length, title_text, ...
                S.fading_param_filename);
            
            % enable load sensor position file and fading parameter file
            set(S.fm,{'enable'},{'on';'on'});
            
        end
        
    end

    function [] = fh_main_crfcn(varargin)
        % Closerequestfcn for figures.
        delete(S.fh) % Delete all figures stored in structure.
    end

    function [] = fh_sub_crfcn(varargin)
        % Closerequestfcn for figures.
        delete(S.fh(2)) % Delete all figures stored in structure.
        
        % #############################
        % ### MUST USE THIS
        % #############################
        S.fh(2) = [];
    end

end
