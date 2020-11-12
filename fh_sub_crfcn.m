function [] = fh_sub_crfcn(varargin)

hObject = varargin{3};
S = guidata(hObject);

% Closerequestfcn for figures.
delete(S.fh(2)) % Delete all figures stored in structure.

% #############################
% ### MUST USE THIS
% #############################
S.fh(2) = 0;

guidata(hObject, S);

end

