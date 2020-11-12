function [] = learn_try_catch(filename)
%
% [usage]
% learn_try_catch('Van_Gogh_Starry_Night.jpg')
% learn_try_catch('Van_Gogh_Starry_Night.jpeg')

imageData = readImage(filename);

end

%%
function imageData = readImage(filename)

try
    imageData = imread(filename);
catch exception
    
    % Is the error because MATLAB could not find the file?
    if ~exist(filename, 'file')
        
        % Check for common typos in the extension.
        [~, ~, extension] = fileparts(filename);
        switch extension
            case '.jpg'
                altFilename = strrep(filename, '.jpg', '.jpeg');
            case '.jpeg'
                altFilename = strrep(filename, '.jpeg', '.jpg');
            case '.tif'
                altFilename = strrep(filename, '.tif', '.tiff');
            case '.tiff'
                altFilename = strrep(filename, '.tiff', '.tif');
            otherwise
                rethrow(exception);
        end
        
        % Try again, with modifed filename.
        try
            imageData = imread(altFilename);
            disp('altfile');
        catch exception2
            % Rethrow original error.
            rethrow(exception)
        end
        
    else
        rethrow(exception)
    end
    
end

end