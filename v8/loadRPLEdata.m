function [data, err] = loadRPLEdata(path)
%% loadRPLEdata.m
% Load in all of the RPLE data contained in the subfolders of the specified
% directory (path). If any of the data has not yet been processed, have the
% user process it using the single spectrum GUI.
    
%% 
    % Let calling program know if there is an error
    err = 0;
    
    % Load in the structure for the folder with all of the data folders
    directory = dir(path);
    
    % Create data structure to store data and metadata for each scan
    chars = char(1:length(directory));
    data = struct('filename',cellstr((chars(1:length(directory)-2))'),...
        'x',cellstr((chars(1:length(directory)-2))'),...
        'y',cellstr((chars(1:length(directory)-2))'),...
        'sy',cellstr((chars(1:length(directory)-2))'),...
        'ACvalue',cellstr((chars(1:length(directory)-2))'),...
        'magvalue',cellstr((chars(1:length(directory)-2))'),...
        'detPol',cellstr((chars(1:length(directory)-2))'),...
        'resOD',cellstr((chars(1:length(directory)-2))'),...
        'HeNeOD',cellstr((chars(1:length(directory)-2))'));
    
    % For loop to load in the spectrum and assign them to the 'data'
    % structure sequentially. Because there may be other files in this folder
    % there will be extra spaces in the structure.
    once = 0;
    for i = 3:length(directory)
        if directory(i).isdir == 1
            subfolder = dir([path '\' directory(i).name]);
            spectrumfound = 0; % Local feedback
            for j = 3:length(subfolder) % Search the subfolder for any files that have 'spectrum' in them
                split = strsplit(subfolder(j).name, ' ');
                isspectrumA = ismember(split, 'spectrum');
                isspectrumB = ismember(split, 'Spectrum');
                if any( [isspectrumA isspectrumB] )
                    k = j;
                    spectrumfound = 1;
                end
            end
            
            if spectrumfound == 0 % If no spectrum file is found
                cprintf('err', ['\nERROR: There is no spectrum file in ' directory(i).name ' .\n']);
                err = 1;
                beep; return
            end
            
            isUNprocessed = ismember(split, '(unprocessed)');
            if any(isUNprocessed)
                if once == 0
                    fprintf(1, '\nFinal processing done on the following scans:\n');
                    once = 1;
                end
                % Data not fully processed, open GUI to let user finish
                for l = 3:length(subfolder)
                    split1a = strsplit(subfolder(l).name, '.');
                    isfigure = ismember(split1a, 'fig');
                    if any(isfigure) % Open the figure
                        open([ path '\' directory(i).name '\' subfolder(l).name ])
                        uiwait()
                        
                        % Go back through and delete the unprocessed files
                        newsubfolder = dir([path '\' directory(i).name]);
                        for m = length(newsubfolder):-1:3
                            split1b = strsplit(newsubfolder(m).name, ' ');
                            isUNprocessed = ismember(split1b, '(unprocessed)');
                            if any(isUNprocessed)
                                delete([ path '\' directory(i).name '\' newsubfolder(m).name ])
                                newsubfolder = dir([path '\' directory(i).name]);
                            end
                        end
                    end
                end
                subfolder = newsubfolder;
                pause(3)
            end
            % Data now fully analyzed, so load it in
            load([ path '\' directory(i).name '\' subfolder(k).name ], 'x', 'y', 'sy');
            data(i-2).filename = directory(i).name; % Loads in the filename (foldername)
            if x(1) < 3000 % Data in wavelength
                data(i-2).x = 299792458./x; % Convert to frequency
            else % Data in frequency
                data(i-2).x = x;
            end
            data(i-2).y = y;
            data(i-2).sy = sy;
        else
            continue,
        end
    end
    
    % This loop goes through in reverse order to delete the parts of the 'data'
    % structure that were left unused.
    for i = length(data):-1:1
        if ischar(data(i).x) == 1
            data(i)=[];
        end
    end
    
end