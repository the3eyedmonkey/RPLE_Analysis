function [tolFound, usertolerance] = ranBefore(path)
%% ranBefore.m
% Determines if the RPLE_Analysis code has been ran before for this set of
% data.

%% 
    % Load in the structure for the folder with all of the data folders.
    directory = dir(path);
    
    % Turn warnings off so matlab does not spit one if the variables below are not found
    warning('off', 'all')
    
    tolFound = 0;
    for i = 3:length(directory)
        split2a = strsplit(directory(i).name, ' ');
        if length(split2a) > 1 %If there are no spaces split2a{end-1}=0 and there is an error, so avoid this
            if strcmp([split2a{end-1} ' ' split2a{end}], 'RPLE data.mat')
                load([ path '\' directory(i).name ], 'usertolerance');
                if exist('usertolerance', 'var')
                    tolFound = 1;
                end
            else
                continue;
            end
        end
    end
    
    % Turn warnings back on again
    warning('on', 'all')
    
end