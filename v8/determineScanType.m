function scanType = determineScanType(data)
%% determineScanType.m
% This program determines whether the scans were AC Stark scans(ref, AC,
% mag, ACmag). Or if they were not AC Stark scans. Additionally
% functionality can be added here by specifying more types of scans, and
% added analysis for those later in the program.

%% 
    % scanType: 0 means reference, 1 means AC, 2 means magnet, 3 means AC
    % magnet, and 4 means other at the moment.
    scanType = zeros(1, length(data));
    for i = 1:length(data)
        AC = data(i).ACvalue;
        mag = data(i).magvalue;
        if AC == 0 && mag == 0
            scanType(i) = 0; % reference
        elseif isnumeric(AC) && (mag == 0 || ~isnumeric(mag))
            scanType(i) = 1; % AC
        elseif isnumeric(mag) && (AC == 0 || ~isnumeric(AC))
            scanType(i) = 2; % magnet
        elseif isnumeric(AC) && isnumeric(mag)
            scanType(i) = 3; % AC magnet
        end
        
        
    end

end