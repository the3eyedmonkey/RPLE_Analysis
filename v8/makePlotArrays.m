function plotarrays = makePlotArrays(data, fitvalues, scanType, numPeaks)
%% makePlotArrays.m
% Create plot arrays structure to inform on one thing about all scans, 
% specifically for AC Stark scans.

%% 
    % Total number of scans
    N = length(fitvalues);
    
    % Create AC plotarrays structure
    chars = char(1:N);
    plotarrays = struct('indices',struct('ref',chars(1),'AC',chars(1),'mag',chars(1),'ACmag',chars(1),...
            'res',chars(1),'HeNe',chars(1)),...
        'ACvalues',struct('AC',chars(1),'ACmag',chars(1)),...
        'magvalues',struct('mag',chars(1),'ACmag',chars(1)),...
        'resODs',chars(1),...
        'HeNeODs',chars(1),...
        'w0s',struct('ref',chars(1),'referr',chars(1),'AC',chars(1),'ACerr',chars(1),...
            'mag',chars(1),'magerr',chars(1),'ACmag',chars(1),'ACmagerr',chars(1),...
            'res',char(1),'reserr',chars(1),'HeNe',chars(1),'HeNeerr',chars(1)),...
        'linewidths',struct('ref',chars(1),'referr',chars(1),'AC',chars(1),'ACerr',chars(1),...
            'mag',chars(1),'magerr',chars(1),'ACmag',chars(1),'ACmagerr',chars(1),...
            'res',char(1),'reserr',chars(1),'HeNe',chars(1),'HeNeerr',chars(1)),...
        'heights',struct('ref',chars(1),'referr',chars(1),'AC',chars(1),'ACerr',chars(1),...
            'mag',chars(1),'magerr',chars(1),'ACmag',chars(1),'ACmagerr',chars(1),...
            'res',char(1),'reserr',chars(1),'HeNe',chars(1),'HeNeerr',chars(1)),...
        'areas',struct('ref',chars(1),'referr',chars(1),'AC',chars(1),'ACerr',chars(1),...
            'mag',chars(1),'magerr',chars(1),'ACmag',chars(1),'ACmagerr',chars(1),...
            'res',char(1),'reserr',chars(1),'HeNe',chars(1),'HeNeerr',chars(1)),...
        'B',struct('ref',chars(1),'referr',chars(1),'AC',chars(1),'ACerr',chars(1),...
            'mag',chars(1),'magerr',chars(1),'ACmag',chars(1),'ACmagerr',chars(1),...
            'res',char(1),'reserr',chars(1),'HeNe',chars(1),'HeNeerr',chars(1)));
    
    % Initialize arrays
    fn1 = fieldnames(plotarrays);
    for i = 1:numel(fn1)
        if i ~= 4 && i ~= 5
            fn2 = fieldnames(plotarrays.(fn1{i}));
            for j = 1:numel(fn2)
                plotarrays.(fn1{i}).(fn2{j}) = zeros(1, 2*N);
            end
        else
            plotarrays.(fn1{i}) = zeros(1, 2*N);
        end
    end
    
    % For loop to run through all scans
    for i = 1:N
        j = i + (N-1); % Index for if 2 Lorentzians
        % Check scan type
        if scanType(i) == 0 % Reference scan
            type = 'ref';
            typeerr = 'referr';
        elseif scanType(i) == 1 % AC scan
            type = 'AC';
            typeerr = 'ACerr';
            plotarrays.ACvalues.AC(i) = data(i).ACvalue;
            if numPeaks(i) == 2 % Fit with 2 Lorentzians
                plotarrays.ACvalues.AC(j) = data(i).ACvalue;
            end
        elseif scanType(i) == 2 % Magnet scan
            type = 'mag';
            typeerr = 'magerr';
            plotarrays.magvalues.mag(i) = data(i).magvalue;
            if numPeaks(i) == 2 % Fit with 2 Lorentzians
                plotarrays.magvalues.mag(j) = data(i).magvalue;
            end
        elseif scanType(i) == 3 % ACmagnet scan
            type = 'ACmag';
            typeerr = 'ACmagerr';
            plotarrays.ACvalues.ACmag(i) = data(i).ACvalue;
            plotarrays.magvalues.ACmag(i) = data(i).magvalue;
            if numPeaks(i) == 2 % Fit with 2 Lorentzians
                plotarrays.ACvalues.ACmag(j) = data(i).ACvalue;
                plotarrays.magvalues.ACmag(j) = data(i).magvalue;
            end
        elseif scanType(i) == 4 % Resonant OD scan
            type = 'res';
            typeerr = 'reserr';
            plotarrays.resODs(i) = data(i).resOD;
            if numPeaks(i) == 2 % Fit with 2 Lorentzians
                plotarrays.resODs(j) = data(i).resOD;
            end
        elseif scanType(i) == 5 % HeNe OD scan
            type = 'HeNe';
            typeerr = 'HeNeerr';
            plotarrays.resODs(i) = data(i).resOD;
            if numPeaks(i) == 2 % Fit with 2 Lorentzians
                plotarrays.HeNeODs(j) = data(i).resOD;
            end
        end
        
        % Assign to arrays
        plotarrays.(fn1{1}).(type)(i) = i; % Scan index
        
        % Fit parameters
        plotarrays.(fn1{6}).(type)(i) = fitvalues(i).w0s(1);
        plotarrays.(fn1{7}).(type)(i) = fitvalues(i).linewidths(1);
        plotarrays.(fn1{8}).(type)(i) = fitvalues(i).heights(1);
        plotarrays.(fn1{9}).(type)(i) = fitvalues(i).areas(1);
        plotarrays.(fn1{10}).(type)(i) = fitvalues(i).B(1);
        
        % Fit parameter confidence intervals (symmetric so upper-lower gives error bar)
        plotarrays.(fn1{6}).(typeerr)(i) = fitvalues(i).w0s(3)-fitvalues(i).w0s(2);
        plotarrays.(fn1{7}).(typeerr)(i) = fitvalues(i).linewidths(3)-fitvalues(i).linewidths(2);
        plotarrays.(fn1{8}).(typeerr)(i) = fitvalues(i).heights(3)-fitvalues(i).heights(2);
        plotarrays.(fn1{9}).(typeerr)(i) = fitvalues(i).areas(3)-fitvalues(i).areas(2);
        plotarrays.(fn1{10}).(typeerr)(i) = fitvalues(i).B(3)-fitvalues(i).B(2);
        
        if numPeaks(i) == 2 % Fit with 2 Lorentzians
            plotarrays.(fn1{1}).(type)(j) = i; % Scan index

            % Fit parameters
            plotarrays.(fn1{6}).(type)(j) = fitvalues(i).w0s(4);
            plotarrays.(fn1{7}).(type)(j) = fitvalues(i).linewidths(4);
            plotarrays.(fn1{8}).(type)(j) = fitvalues(i).heights(4);
            plotarrays.(fn1{9}).(type)(j) = fitvalues(i).areas(4);
            plotarrays.(fn1{10}).(type)(j) = plotarrays.(fn1{8}).(type)(i); % Only one background parameter

            % Fit parameter confidence intervals (symmetric so upper-lower gives error bar)
            plotarrays.(fn1{6}).(typeerr)(j) = fitvalues(i).w0s(6)-fitvalues(i).w0s(5);
            plotarrays.(fn1{7}).(typeerr)(j) = fitvalues(i).linewidths(6)-fitvalues(i).linewidths(5);
            plotarrays.(fn1{8}).(typeerr)(j) = fitvalues(i).heights(6)-fitvalues(i).heights(5);
            plotarrays.(fn1{9}).(typeerr)(j) = fitvalues(i).areas(6)-fitvalues(i).areas(5);
            plotarrays.(fn1{10}).(typeerr)(j) = plotarrays.(fn1{8}).(typeerr)(i); % Only one background parameter
        end
    end
    
    % Go backwards through arrays deleting unassigned
    for i = 2*N:-1:1
        for j = 1:numel(fn1)
            if j ~= 4 && j ~= 5
                fn2 = fieldnames(plotarrays.(fn1{j}));
                for k = 1:numel(fn2)
                    if plotarrays.(fn1{j}).(fn2{k})(i) == 0
                        plotarrays.(fn1{j}).(fn2{k})(i) = [];
                    end
                end
            else
                if plotarrays.(fn1{j})(i) == 0
                    plotarrays.(fn1{j})(i) = [];
                end
            end
        end
    end
    
    % Some assertions to assure that the structure was made correctly
    % Check that the reference arrays are correct length
    if length(plotarrays.indices.ref) ~= length(plotarrays.w0s.ref)
        cprintf('err', '\nERROR: The length of reference plot arrays do not match.\n');
        beep; return
    end
    
    % Check that the AC arrays are correct length
    if length(plotarrays.indices.AC) ~= length(plotarrays.w0s.AC)
        cprintf('err', '\nERROR: The length of AC plot arrays do not match.\n');
        beep; return
    end
    
    % Check that the magnet arrays are the correct length
    if length(plotarrays.indices.mag) ~= length(plotarrays.w0s.mag)
        cprintf('err', '\nERROR: The length of magnet plot arrays do not match.\n');
        beep; return
    end
    
    % Check that the ACmagnet arrays are the correct length
    if length(plotarrays.indices.ACmag) ~= length(plotarrays.w0s.ACmag)
        cprintf('err', '\nERROR: The length of ACmagnet plot arrays do not match.\n');
        beep; return
    end
    
end