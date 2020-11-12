function [] = fh_main_crfcn(varargin)

hObject = varargin{1};
S = guidata(hObject);

% ########### final solution:
% ########### (1) use try, catch
% ########### or
% ########### (2) check ishghandle(S.fh(2)) && S.fh(2) ~= 0 (current implementation)

%         % Closerequestfcn for figures.
%         % ######## use try, catch to circumvent error, below comments
%         try
%             delete(S.fh); % Delete all figures stored in structure.
%         catch
%             ;
%         end

% #############################################################################################
% ##### when sub figure was not created, closing main figure make error:
% ##### [matlab error message]
%
% sub fig valid
% Error using delete
% Root object may not be deleted
% Error in gui_sim_tdoa/fh_main_crfcn (line 1260)
%          delete(S.fh(2));
% Error while evaluating figure CloseRequestFcn
%
% ##### after sub figure created was closed, closing main figure DONT make error
% #####
% ##### when sub figure was not created even, why ishghandle(S.fh(2)) return with "true"?
% ##### [answer?] 0 fig handle is root object, so ishghandle NOT work for root object
% #############################################################################################

% Closerequestfcn for figures.
delete(S.fh(1)); % Delete main figure stored in structure.

if ishghandle(S.fh(2)) && S.fh(2) ~= 0
    %             disp('sub fig valid');
    delete(S.fh(2)); % delete sub figure
end

% ############ S.fh is no more valid
% ############ DONT need update gui shared data, this program is about to exiting
% guidata(hObject, S);

end

