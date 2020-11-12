function [answer] = learn_inputdlg

prompt = {'Enter matrix size:','Enter colormap name:'};
dlg_title = 'Input';

% (1) If num_lines is a scalar, it applies to all prompts.
% (2) If num_lines is a column vector, each element specifies the number of lines of input for a prompt.
% (3) If num_lines is an array, it must be size m-by-2, where m is the number of prompts on the dialog box. 
%     Each row refers to a prompt. 
%     The first column specifies the number of lines of input for a prompt. 
%     The second column specifies the width of the field in characters.
num_lines = [2; 5];
def = {'20','hsv'}; % default
answer = inputdlg(prompt,dlg_title,num_lines,def);

end