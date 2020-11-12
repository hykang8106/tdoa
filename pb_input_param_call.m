function [] = pb_input_param_call(varargin)

hObject = varargin{1};
S = guidata(hObject);
% ###################################
% ### in here, S is good as expected
% ###################################
S;

% disp('here pb input');
% ###### function FAILED!!!!!!!
%         sub_pb_input_param_call(S);

% Callback for the pushbutton.
sel = findobj(get(S.bg, 'selectedobject'));  % See BUG note in GUI_8
S.fh;
S.SEL_old;
S.SEL = find(S.rd == sel);  % Store current radiobutton.
if (S.fh(2) ~= 0) && (S.SEL ~= S.SEL_old)
    delete(S.fh(2));
    S.fh(2) = 0;
end
S.SEL_old = S.SEL;
S.SEL;

switch sel
    
    % ##########################
    % ##### sensor fixed
    % ##########################
    case S.rd(1)
        
        if S.fh(2) == 0  % We haven't been here before.
            
            % ####### failed to make this code function
            % ####### error message: can't find callback function
            % ####### why?! rewrite code
            %                     gui_populate_sensor_fixed_uicontrol(S);
            
            [S] = gui_component_layout_sensor_fixed(S);
            
%             S.fh(2) = figure('units','pixels',...
%                 'position',[S.second_fig_pos_left 500 350 430],... % [1070 500 350 430]
%                 'name','input for sensor fixed',...
%                 'menubar','none',...
%                 'numbertitle','off',...
%                 'resize','off',...
%                 'closerequestfcn',{@fh_sub_crfcn, S.fh(1)});
%             
%             S.tx(1) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[70 385 210 25],...
%                 'string','input for sensor fixed',...
%                 'fontweight','bold',...
%                 'fontsize',10, ...
%                 'foregroundcolor', 'b', ...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.ed(1) = uicontrol('units','pixels',...
%                 'style','edit',...
%                 'unit','pix',...
%                 'string', '100', ...
%                 'position',[170 210 150 25]);
%             S.tx(2) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[30 205 110 25],...
%                 'string','trial length',...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.ed(2) = uicontrol('units','pixels',...
%                 'style','edit',...
%                 'unit','pix',...
%                 'string', '1.5', ...
%                 'position',[170 245 150 25]);
%             S.tx(3) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[30 240 110 25],...
%                 'string','radius ratio',...
%                 'tooltipstring', 'distance ratio between target and farthest sensor from [0,0]', ...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.ed(3) = uicontrol('units','pixels',...
%                 'style','edit',...
%                 'unit','pix',...
%                 'string', '100', ...
%                 'position',[170 280 150 25]);
%             S.tx(4) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[30 275 110 25],...
%                 'string','y of target position',...
%                 'tooltipstring', 'unit = meter', ...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.ed(4) = uicontrol('units','pixels',...
%                 'style','edit',...
%                 'unit','pix',...
%                 'string', '100', ...
%                 'position',[170 315 150 25]);
%             S.tx(5) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[30 310 110 25],...
%                 'string','x of target position',...
%                 'tooltipstring', 'unit = meter', ...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.ed(5) = uicontrol('units','pixels',...
%                 'style','edit',...
%                 'unit','pix',...
%                 'string', '5', ...
%                 'position',[170 350 150 25]);
%             S.tx(6) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[30 345 110 25],...
%                 'string','SNR in dB',...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.pb(2) = uicontrol('style','push',...
%                 'unit','pix',...
%                 'position',[100 20 150 30],...
%                 'fontweight','bold',...
%                 'string','run sensor fixed',...
%                 'callback',{@pb_run_sensor_fixed_call, S.fh(1)});
            
        end
        
        % ##########################
        % ##### target fixed
        % ##########################
    case S.rd(2)
        
        if S.fh(2) == 0  % We haven't been here before.
            
            [S] = gui_component_layout_target_fixed(S);
            
%             S.fh(2) = figure('units','pixels',...
%                 'position',[S.second_fig_pos_left 298 417 552],... % [520 298 417 552]
%                 'name','input for target fixed',...
%                 'menubar','none',...
%                 'numbertitle','off',...
%                 'resize','off',...
%                 'closerequestfcn',{@fh_sub_crfcn, S.fh(1)});
%             
%             S.tx(1) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[150 502 150 15],...
%                 'string','input for target fixed',...
%                 'fontweight','bold',...
%                 'fontsize',10, ...
%                 'foregroundcolor', 'b', ...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.lb(1) = uicontrol('units','pixels',...
%                 'style','listbox',...
%                 'unit','pix',...
%                 'string', {'3', '4', '5', '6', '7'}, ...
%                 'value', 3, ... % default: index = 3, sensor length = 5
%                 'position',[176 377 169 75]);
%             S.tx(2) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[69 410 91 15],...
%                 'string','sensor length',...
%                 'HorizontalAlignment', 'left', ...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.ed(1) = uicontrol('units','pixels',...
%                 'style','edit',...
%                 'unit','pix',...
%                 'string', '5', ...
%                 'position',[176 329 52 23]);
%             S.tx(3) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[68 333 77 15],...
%                 'string','SNR in dB',...
%                 'HorizontalAlignment', 'left', ...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.cb(1) = uicontrol('units','pixels',...
%                 'style','checkbox',...
%                 'unit','pix',...
%                 'string', 'sensor 1 is always nearest to target', ...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'), ...
%                 'value', 0, ...
%                 'position',[50 279 244 23]);
%             S.cb(2) = uicontrol('units','pixels',...
%                 'style','checkbox',...
%                 'unit','pix',...
%                 'string', 'use only least square estimator', ...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'), ...
%                 'value', 0, ...
%                 'position',[50 229 217 23]);
%             S.cb(3) = uicontrol('units','pixels',...
%                 'style','checkbox',...
%                 'unit','pix',...
%                 'string', 'plot sensor position & hyperbola curve', ...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'), ...
%                 'value', 1, ...
%                 'position',[50 179 257 23]);
%             S.cb(4) = uicontrol('units','pixels',...
%                 'style','checkbox',...
%                 'unit','pix',...
%                 'string', 'plot target signal & correlation', ...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'), ...
%                 'value', 0, ...
%                 'position',[50 129 208 23]);
%             
%             S.pb(2) = uicontrol('style','push',...
%                 'unit','pix',...
%                 'position',[137 39 114 23],...
%                 'fontweight','bold',...
%                 'string','run target fixed',...
%                 'callback',{@pb_run_target_fixed_call, S.fh(1)});
            
        end
        
        % ##########################
        % ##### batch sensor fixed
        % ##########################
    case S.rd(3)
        if S.fh(2) == 0  % We haven't been here before.
            
            [S] = gui_component_layout_batch_sensor_fixed(S);
            
%             S.fh(2) = figure('units','pixels',...
%                 'position',[S.second_fig_pos_left 500 350 430],... % [1070 500 350 430]
%                 'name','input for batch sensor fixed',...
%                 'menubar','none',...
%                 'numbertitle','off',...
%                 'resize','off',...
%                 'closerequestfcn',{@fh_sub_crfcn, S.fh(1)});
%             
%             S.tx(1) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[70 385 210 25],...
%                 'string','input for batch sensor fixed',...
%                 'fontweight','bold',...
%                 'fontsize',10, ...
%                 'foregroundcolor', 'b', ...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.ed(1) = uicontrol('units','pixels',...
%                 'style','edit',...
%                 'unit','pix',...
%                 'string', '100', ...
%                 'position',[170 210 150 25]);
%             S.tx(2) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[30 205 110 25],...
%                 'string','trial length',...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.ed(2) = uicontrol('units','pixels',...
%                 'style','edit',...
%                 'unit','pix',...
%                 'string', '1.5', ...
%                 'position',[170 245 150 25]);
%             S.tx(3) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[30 240 110 25],...
%                 'string','radius ratio',...
%                 'tooltipstring', 'ratio between target distance and farthest sensor distance from [0,0]', ...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.ed(4) = uicontrol('units','pixels',...
%                 'style','edit',...
%                 'unit','pix',...
%                 'string', '100', ...
%                 'position',[170 315 150 25]);
%             S.tx(5) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[30 310 110 25],...
%                 'string','target position step',...
%                 'tooltipstring', 'same step is applied to x and y, unit = meter', ...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.ed(5) = uicontrol('units','pixels',...
%                 'style','edit',...
%                 'unit','pix',...
%                 'string', '5', ...
%                 'position',[170 350 150 25]);
%             S.tx(6) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[30 345 110 25],...
%                 'string','SNR in dB',...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.pb(2) = uicontrol('style','push',...
%                 'unit','pix',...
%                 'position',[100 20 150 30],...
%                 'fontweight','bold',...
%                 'string','run batch sensor fixed',...
%                 'callback',{@pb_run_batch_sensor_fixed_call, S.fh(1)});
            
        end
        
        % ##########################
        % ### batch target fixed
        % ##########################
    case S.rd(4)
        
        if S.fh(2) == 0  % We haven't been here before.
            
            [S] = gui_component_layout_batch_target_fixed(S);
            
%             S.fh(2) = figure('units','pixels',...
%                 'position',[S.second_fig_pos_left 298 417 552],... % [520 298 417 552]
%                 'name','input for target fixed',...
%                 'menubar','none',...
%                 'numbertitle','off',...
%                 'resize','off',...
%                 'closerequestfcn',{@fh_sub_crfcn, S.fh(1)});
%             
%             S.tx(1) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[90 501 185 19],...
%                 'string','input for batch target fixed',...
%                 'fontweight','bold',...
%                 'fontsize',10, ...
%                 'foregroundcolor', 'b', ...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.ed(1) = uicontrol('units','pixels',...
%                 'style','edit',...
%                 'unit','pix',...
%                 'string', '3:7', ...
%                 'position',[177   405    72    23]);
%             S.tx(2) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[32 410 128 15],...
%                 'string','sensor length (vector)',...
%                 'HorizontalAlignment', 'left', ...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.ed(2) = uicontrol('units','pixels',...
%                 'style','edit',...
%                 'unit','pix',...
%                 'string', '-10:5:15', ...
%                 'position',[176 329 72 23]);
%             S.tx(3) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[32 333 112 15],...
%                 'string','SNR in dB (vector)',...
%                 'HorizontalAlignment', 'left', ...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.ed(3) = uicontrol('units','pixels',...
%                 'style','edit',...
%                 'unit','pix',...
%                 'string', '100', ...
%                 'position',[176   254    72    23]);
%             S.tx(4) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[32   258   111    15],...
%                 'string','trial length',...
%                 'HorizontalAlignment', 'left', ...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.pb(2) = uicontrol('style','push',...
%                 'unit','pix',...
%                 'position',[115 39 175 23],...
%                 'fontweight','bold',...
%                 'string','run batch target fixed',...
%                 'callback',{@pb_run_batch_target_fixed_call, S.fh(1)});
            
        end
        
        % ##########################
        % ### batch bw corr snr
        % ##########################
    case S.rd(5)
        
        if S.fh(2) == 0  % We haven't been here before.
            
            [S] = gui_component_layout_batch_bw_corr_snr(S);
            
%             S.fh(2) = figure('units','pixels',...
%                 'position',[S.second_fig_pos_left   290   361   510],... % [520   290   361   510]
%                 'name','input for batch bw corr snr',...
%                 'menubar','none',...
%                 'numbertitle','off',...
%                 'resize','off',...
%                 'closerequestfcn',{@fh_sub_crfcn, S.fh(1)});
%             
%             S.tx(1) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[100   477   151    15],...
%                 'string','input for batch bw corr snr',...
%                 'fontweight','bold',...
%                 'fontsize',10, ...
%                 'foregroundcolor', 'b', ...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.ed(1) = uicontrol('units','pixels',...
%                 'style','edit',...
%                 'unit','pix',...
%                 'string', '100, 100', ...
%                 'position',[200   437   101    23]);
%             S.tx(2) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[50   441   120    15],...
%                 'string','target position (x, y)',...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.ed(2) = uicontrol('units','pixels',...
%                 'style','edit',...
%                 'unit','pix',...
%                 'string', '1.5', ...
%                 'position',[200   387   101    23]);
%             S.tx(3) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[50   391   101    15],...
%                 'string','radius ratio',...
%                 'tooltipstring', 'distance ratio between target and farthest sensor from [0,0]', ...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.ed(3) = uicontrol('units','pixels',...
%                 'style','edit',...
%                 'unit','pix',...
%                 'string', '100', ...
%                 'position',[200   337   101    23]);
%             S.tx(4) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[50   341   101    15],...
%                 'string','trial length',...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.lb(1) = uicontrol('units','pixels',...
%                 'style','listbox',...
%                 'unit','pix',...
%                 'string', {'6', '15', '25'}, ...
%                 'value', 2, ...
%                 'position',[200   242   101    68]);
%             S.tx(5) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[50   273   101    15],...
%                 'string','NDLRB',...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.ed(4) = uicontrol('units','pixels',...
%                 'style','edit',...
%                 'unit','pix',...
%                 'string', '-5, 0, 5', ...
%                 'position',[200   187   101    23]);
%             S.tx(6) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[50   191   101    15],...
%                 'string','SNR in dB (vector)',...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.ed(5) = uicontrol('units','pixels',...
%                 'style','edit',...
%                 'unit','pix',...
%                 'string', '1:15', ...
%                 'position',[200   137   101    23]);
%             S.tx(7) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[50   141   101    15],...
%                 'string','NPRSRB (vector)',...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.ed(6) = uicontrol('units','pixels',...
%                 'style','edit',...
%                 'unit','pix',...
%                 'string', '1, 2', ...
%                 'position',[200    87   101    23]);
%             S.tx(8) = uicontrol('units','pixels',...
%                 'style','text',...
%                 'unit','pix',...
%                 'position',[50    91   101    15],...
%                 'string','nsubframe (vector)',...
%                 'fontweight','bold',...
%                 'backgroundcolor',get(S.fh(2),'color'));
%             
%             S.pb(2) = uicontrol('style','push',...
%                 'unit','pix',...
%                 'position',[100    37   151    23],...
%                 'fontweight','bold',...
%                 'string','run batch bw corr snr',...
%                 'callback',{@pb_run_batch_bw_corr_snr_call, S.fh(1)});
            
        end
        
        % ##########################
        % ### plot target fixed
        % ##########################
    case S.rd(6)
        
        filterspec = 'tdoa_result_target_fixed_*.mat';
        [filename, pathname, filterindex] = uigetfile(filterspec);
        if ~filename
            fprintf(2, '######## file selection canceled\n');
            return;
        else
            fprintf('filename = %s\n', filename);
        end
        
        try
            V = load(filename, 'sensor_length');
            prompt = {sprintf('sensor length for errorbar plot (one of %s):', num2str(V.sensor_length))};
            dlg_title = 'Input';
            num_lines = [1];
            def = {'5'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            sensor_number_torrieri = str2num(answer{1});
        catch
            disp('Unable to Load. Check Name and Try Again.');
        end
        
        sub_plot_tdoa_result_target_fixed(filename, sensor_number_torrieri);
        
        % ##########################
        % ### plot sensor fixed
        % ##########################
    case S.rd(7)
        
        filterspec = 'tdoa_result_sensor_fixed_*.mat';
        [filename, pathname, filterindex] = uigetfile(filterspec);
        if ~filename
            fprintf(2, '######## file selection canceled\n');
            return;
        else
            fprintf('filename = %s\n', filename);
        end
        
        sub_plot_tdoa_result_sensor_fixed(filename);
        
        % ##########################
        % ### plot bw corr snr
        % ##########################
    case S.rd(8)
        
        filterspec = 'tdoa_result_bw_corr_snr_*.mat';
        [filename, pathname, filterindex] = uigetfile(filterspec);
        if ~filename
            fprintf(2, '######## file selection canceled\n');
            return;
        end
        
        try
            % ### reminder ###
            % ### in batch_simulate_tdoa_bw_corr_snr.m
            % save(filename, 'sensor_position', 'target_position', 'trial_length', 'rician_param', ...
            %     'ndlrb', 'snr_vec', 'nprsrb_vec', 'nsubframe_vec', 'position_error_array', 'target_position_limit', ...
            %     'error_cell_array', 'radius_ratio');
            
            % old result mat file may not have 'rician_param' variable
            V = load(filename, 'rician_param', 'snr_vec', 'nsubframe_vec');
            if ~isempty(V.rician_param)
                prompt = {sprintf('index for SNR = [%s] dB:', num2str(V.snr_vec)), ...
                    sprintf('index for nsubframe = [%s]:', num2str(V.nsubframe_vec))};
                dlg_title = 'Input';
                num_lines = 1;
                def = {'1', '1'};
                answer = inputdlg(prompt,dlg_title,num_lines,def);
                snr_idx = str2num(answer{1});
                nsubframe_idx = str2num(answer{2});
            else
                % ### below 2 line is dummy for function input in case of non-fading environment
                snr_idx = 1;
                nsubframe_idx = 1;
            end
            
            sub_plot_tdoa_result_bw_corr_snr(filename, snr_idx, nsubframe_idx);
            
        catch
            disp('Unable to Load. Check Name and Try Again.');
        end
        
    otherwise
        % Very unlikely I think.
        
end

guidata(hObject, S);

end

