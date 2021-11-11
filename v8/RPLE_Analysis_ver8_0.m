function RPLE_Analysis_ver8_0(varargin)
%% RPLE_Analysis_ver8_0.m
% By Tristan Wilkinson
% 
% This program is used for loading and analyzing a bunch of RPLE spectra.
% Specfically comparing reference and AC scans to determine properties of
% the AC Stark shift and Dynamic Nuclear Polarization (DNP).
%  
% Two main structures are created. The fitvalues structure contains all of
% the fitting parameters for each individual scan, and can be parsed by
% using the scan indices. The plotarrays structure contains the same
% information, but instead is organized to contain a single value for all
% scans.
% 
% ******* ver8.0 (11/5/21) **********
% UPDATES:
% 
% This version will be a complete overhaul to allow for generalized
% functionality with the TAW_RPLE_MultipleScans_v2.vi.
% 
% I am also going to attempt to modulaized the code and figure out how to
% use Git to my advantage.
% 
% 
% PROGRAM OVERVIEW:  (goes roughly section by section)
% A spectrum file is looked for in all of the folders in the path
% directory. Each spectrum and parameter settings is loaded into the data
% structure.
% 
% There are 3 overall possibilities for type of scans:
% AC Stark scans: AC Stark laser and/or magnet in use
% Excitation laser power scans: Resonant OD and/or HeNe OD in use
% Other: Any other combinations
% 
% Each scan is fit with a single Lorentzian. If the goodness of fit is <
% tolerance, then the data is fit with a sum of two Lorentzians. If the fit
% is still bad, then the user fits the data manually. The user can also
% choose to visually check each fit and confirm the fit is good. All of the
% fit parameters are placed in the fitvalues structure.
% 
% The structure plotarrays is created, which informs on one thing about all
% scans, whereas the fitvalues structure informs on all things about one
% scan.
% 
% Plots are created for experimental parameters vs. fit parameters. As of
% now I only have built in functionality for AC Stark and excitation
% scans.
% 
% Figures, data, and log file are saved.
% ***********************************

close all
clc
%% Input parameters
    in_params = inputParser;
    in_params.CaseSensitive = false;
    in_params.addParameter('path','', @ischar);
    in_params.parse(varargin{:});
    
    path = in_params.Results.path;
%% Loading in the data
    % Check the path
    if isempty(path)
        % Ask the user to select the folder of interest
        path = uigetdir(\\ecas.wvu.edu\squol\AC Stark Effect);
        %path = '\\ecas.wvu.edu\squol\Data\RPLE_Testing_Fall_2021\MultiScanTesting\RaMultiScanTesting2';
        if path == 0 % User pressed cancel
            cprintf('err', '\nCANCELLED: Folder path selection cancelled.\n');
            return
        end
    end
    
    % Extract and print the name of the folder the user has selected
    split = strsplit(path, '\');
    folder = split(end);
    % Use diary to record the command line output to a text file
    dfile  = [path '\' folder{1} ' RPLE analysis log.txt'];
    if exist(dfile, 'file')
        delete(dfile)
    end
    diary(dfile)
    fprintf(1, ['\nRPLE_Analysis: ', folder{1}, '\n']);
    
    % Load the data
    [data, err] = loadRPLEdata(path);
    if err
        return
    end
    
    % Define total number of scans as N
    N = length(data);
    
    fprintf(1, ['\n' num2str(N) ' RPLE spectra:']);
    
    % Extract the parameter values for each scan
    data = extractParams(data);
    
%% Determine type of scan
    scanType = determineScanType(data);
    
    % Determine the number of types of scans
    % AC Stark scan type
    Nref = sum(ismember(scanType, 0));
    NAC = sum(ismember(scanType, 1));
    Nmag = sum(ismember(scanType, 2));
    NACmag = sum(ismember(scanType, 3));
    % Excitation laser power scan type
    Nres = sum(ismember(scanType, 4));
    NHeNe = sum(ismember(scanType, 5));
    NresHeNe = sum(ismember(scanType, 6));
    % All other scan types
    Nother = sum(ismember(scanType, 7));
    
    % Print number of scans, and determine type of scans
    % First some assertions to check everything makes sense
    if (Nref + NAC + Nmag + NACmag) ~= 0 && (Nres + NHeNe + NresHeNe) ~= 0 && Nother ~= 0
        cprintf('err', '\nERROR: All three types of scans found in the same folder...yikes.\n');
        beep; return
    elseif (Nref + NAC + Nmag + NACmag) ~= 0 && (Nres + NHeNe + NresHeNe) ~= 0
        cprintf('err', '\nERROR: AC Stark and excitation laser power scans found in the same folder.\n');
        beep; return
    elseif (Nref + NAC + Nmag + NACmag) ~= 0 && Nother ~= 0
        cprintf('err', '\nERROR: AC Stark and "other" scans found in the same folder.\n');
        beep; return
    elseif (Nres + NHeNe + NresHeNe) ~= 0 && Nother ~= 0
        cprintf('err', '\nERROR: Excitation laser power and "other" scans found in the same folder.\n');
        beep; return
    elseif (Nres + NHeNe) ~= 0 && NresHeNe ~= 0
        cprintf('err', '\nERROR: Excitation laser power scans do not makes sense.\n');
        beep; return
    elseif (Nref + NAC + Nmag + NACmag) ~= 0
        ACStarkScans = 1;
        excitationLaserPowerScans = 0;
        otherScans = 0;
        fprintf(1, '\tAC Stark scans\n');
        fprintf(1, ['\t\t\t' num2str(Nref) ' reference\n']);
        fprintf(1, ['\t\t\t' num2str(NAC) ' AC\n']);
        fprintf(1, ['\t\t\t' num2str(Nmag) ' magnet\n']);
        fprintf(1, ['\t\t\t' num2str(NACmag) ' AC magnet\n']);
    elseif (Nres + NHeNe + NresHeNe) ~= 0
        ACStarkScans = 0;
        excitationLaserPowerScans = 1;
        otherScans = 0;
        fprintf(1, '\texcitation laser power scans\n');
        fprintf(1, ['\t\t\t' num2str(Nres) ' resonant OD\n']);
        fprintf(1, ['\t\t\t' num2str(NHeNe) ' HeNe OD\n']);
        fprintf(1, ['\t\t\t' num2str(NresHeNe) ' resonant HeNe OD\n']);
    elseif Nother ~= 0
        ACStarkScans = 0;
        excitationLaserPowerScans = 0;
        otherScans = 1;
        fprintf(1, '\tother scans\n');
        fprintf(1, ['\t\t\t' num2str(Nother) ' total\n']);
    end
    
%% Fit with Lorentzian(s)
    % Tolerance for goodness of fit
    tolerance = 0.93;
    
    % Array for feedback on number of peaks in each spectrum
    numPeaks = zeros(1, N);
    
    % Create a structure called 'fitvalues' for all of the fit parameters
    chars = char(1:N);
    fitvalues = struct('filename',cellstr((chars(1:N))'),...
            'w0s',cellstr((chars(1:N))'),...
            'linewidths',cellstr((chars(1:N))'),...
            'heights',cellstr((chars(1:N))'),...
            'areas',cellstr((chars(1:N))'),...
            'B',cellstr((chars(1:N))'));
    
    fprintf(1, '\nFitting progress:\n');
    
    % See if the user wants to visually check each fit for quality
    answer0 = questdlg('Would you like to manually check each fit?',...
        'Manual fitting', 'No');
    
    switch answer0
        case 'No' % Only check if a bad fit
            manualfit = 0;
        case 'Cancel'
            cprintf('err', '\nCANCELLED: Manual fitting selection cancelled.\n');
            return
        case '' % User closed the dialog box
            cprintf('err', '\nCANCELLED: Manual fitting selection cancelled.\n');
            return
        case 'Yes' % Check each fit manually
            manualfit = 1;
            fprintf(1, '\t\t\t\t\tFits checked manually\n');
    end
    
    % For loop to fit all of the spectra
    anymanual = 0;
    for i = 1:N
        
        % Data for fit
        xData = data(i).x;
        yData = data(i).y;
        yErr = data(i).sy;
        
        % In case the error bars are messed up
        if any(yErr == 0)
            yErr = ones(1, length(yErr));
        end
        
        % Fit with a single Lorentzian
        [f, gof, guess] = fitLorentzian(xData, yData, yErr);
        
        % Fit R and L detection with 1 peak only if AC Stark data
        if ACStarkScans && ((data(i).detPol == 'R') || (data(i).detPol == 'L'))
            tol = 0.6;
        else
            tol = tolerance;
        end
        
        % Evaluate rsquared to see if need to fit with 2 Lorentzians
        badfit = 0;
        if gof.rsquare <= tol
            
            % Fit the data with a sum of two Lorentzians
            [f, gof, guess] = fit2Lorentzians(xData, yData, yErr);
            
            % Record if fit is still bad
            if gof.rsquare <= tol
                badfit = 1;
            end
            
            % Set the number of peaks
            numPeaks(i) = 2;
        else
            numPeaks(i) = 1;
        end
        
        % Plot the raw data with fits
        figure('Name', data(i).filename, 'WindowStyle', 'docked', 'numbertitle', 'off');
        hold on
        h0 = errorbar(xData, yData, yErr, 'b.', 'Capsize', 0.1);
        h1 = plot(f);
        
        legend([h0,h1], 'data', 'fit')
        
        xlabel('Frequency (GHz)');
        ylabel('Intensity (arb. units)');
        
        % If fit is bad or user wants to check each one, do manual fitting
        if badfit || manualfit
            anymanual = 1;
            [numPeaks(i), ft, cancelled] = manualFitting(xData, yData, yErr, guess, numPeaks(i));
            if cancelled
                return
            elseif ft
                fprintf(1, ['\t\t\tScan ' data(i).filename ' fit manually\n']);
            end
        end
        
        % Assign parameters to fitvalues structure
        fitvalues(i) = assignFitValues(data(i), fitvalues(i), f, numPeaks(i));
    end
    
    if ~anymanual
        fprintf(1, '\t\t\tAll scans fit automatically\n');
    end
    fprintf(1, ['\t\t\tTolerance: ' num2str(tolerance) '\n']);
    
%% Plots
    % Plot all spectra
    plot_allSpectra(data)
    
    % Create plotarrays structure to contain 1 thing about all scans
    plotarrays = makePlotArrays(data, fitvalues, scanType, numPeaks);
    
    if ACStarkScans % Analysis specific to AC Stark scans
        
        % Mean of reference scans
        refmean = mean(plotarrays.w0s.ref);
        % Print the refmean to command line
        fprintf(1, ['\nThe mean of reference scans is ' num2str(refmean) ' GHz\n']);
        
        % Plot w0s vs. Scan Index
        plot_w0svScanIndex(plotarrays)
        
        % AC plots
        if NAC ~= 0
            plots_AC(plotarrays)
        end
        
        % Magnet plots
        if Nmag ~= 0
            plots_magnet(plotarrays)
        end
        
        % ACmagnet plots
        if NACmag ~= 0
            plots_ACmagnet(plotarrays)
        end
        
    elseif excitationLaserPowerScans % Analysis specific to excitation laser power scans
        
        % Make the plots
        if Nres ~= 0 && NHeNe == 0 % Resonant OD plots
            plots_res(plotarrays)
        elseif NHeNe ~= 0 && Nres == 0 % HeNe OD plots
            plots_HeNe(plotarrays)
        else % Resonant and HeNe plots
            plots_resHeNe(plotarrays)
        end
        
    elseif otherScans
        % Generalized analysis for non-AC scans goes here
        plotarrays = 'I have no built in functionality for this case...sorry :(';
    end
    
%% Save figure window and variables
    % Array with the figure handles
    h = findobj('type', 'figure');
    Nplots = length(h);
    h = 1:Nplots;
    
    % Save the figure
    savefig(h, [ path '\' folder{1} ' plots' '.fig']);
    
    % Save variables from the workspace
    if ACStarkScans
        save([ path '\' folder{1}, ' RPLE data', '.mat'], 'data', 'fitvalues', 'plotarrays',...
            'Nref', 'NAC', 'Nmag', 'NACmag', 'scanType', 'numPeaks', 'tolerance', 'refmean');
    else
        save([ path '\' folder{1}, ' RPLE data', '.mat'], 'data', 'fitvalues', 'plotarrays',...
            'Nother', 'scanType', 'numPeaks', 'tolerance');
    end
    
    fprintf(1, '\nData and figures saved!\n');
    
    % Output the log file
    diary off
    
end