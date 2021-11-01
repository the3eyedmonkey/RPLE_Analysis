function data = extractParams(data)
%% extractParams().m
% This program parses the file names, extracts the parameter values, and
% places them into the data structure.

%% 
    % Run through all scans
    for i = 1:length(data)
        % Split the folder name by spaces
        split = strsplit(data(i).filename, ' ');
        % Iterate through the split
        for j = 1:length(split) % AC Stark power
            if strcmp(split(j), 'AC')
                data(i).ACvalue = str2double(cell2mat(split(j+1)));
            elseif strcmp(split(j), 'magnet') % Magnetic field
                data(i).magvalue = str2double(cell2mat(split(j+1)));
            elseif strcmp(split(j), 'resOD') % Resonant OD
                data(i).resOD = str2double(cell2mat(split(j+1)));
            elseif strcmp(split(j), 'HeNeOD') % HeNe OD
                data(i).HeNeOD = str2double(cell2mat(split(j+1)));
            end
        end
        
        % Detection polarization.
        X = any(ismember(split,'X'));
        Y = any(ismember(split,'Y'));
        D = any(ismember(split,'D'));
        A = any(ismember(split,'A'));
        R = any(ismember(split,'R'));
        L = any(ismember(split,'L'));
        if X
            data(i).detPol = 'X';
        elseif Y
            data(i).detPol = 'Y';
        elseif D
            data(i).detPol = 'D';
        elseif A
            data(i).detPol = 'A';
        elseif R
            data(i).detPol = 'R';
        elseif L
            data(i).detPol = 'L';
        end
    end
    
end