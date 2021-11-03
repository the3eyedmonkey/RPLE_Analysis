function RPLE_Analysis_ver8_0(varargin)
%% RPLE_Analysis_ver8_0.m
% By Tristan Wilkinson
% 
% This program is used for loading and analyzing a bunch of RPLE spectra.
% Specfically comparing reference and AC scans to determine properties of
% the AC Stark shift and Dynamic Nuclear Polarization (DNP).
%  
% Two main structures are created. The fitvalues structure contains all of the
% fitting parameters for each individual scan, and can be parsed by using
% the scan indices. The plotarrays structure contains the same information,
% but instead is organized to contain a single value for all scans.
% 
% ******* ver8.0 (11/2/21) **********
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
% The program asks the user via uigetdir() to specify the
% path to the folder DIRECTLY ABOVE the folders with individual scan data.
% You can change the default folder this search starts in on line 98.
% 
% The program loads all of the spectra into the "data" structure. It does
% this by searching all of the subfolders for a data file with the string
% 'spectrum' in the name.
% 
% The program identifies each scan as either reference or AC by reading the
% name of the folders.
% 
% The program determines if the it has been run before for this set of data
% by searching the directory for the RPLE Analysis Data.m data file.
% 
% The user is asked to enter a value for the tolerance the program will use
% to determine if fits are good enough. The user is also asked if this is
% important data that they would like to fit manually. If the user says no
% the program runs as normal.
%
% Each spectrum is fit with a single Lorentzian, then the gof.rsquared
% value is evaluated for AC scans to determine if fitting with a sum of two
% Lorentzians is warrented. This is dependent on detection polarization; X
% and Y are fit with a double, R and L are fit witht a single.
%
% If the user has indicated that this is important data, then the program
% asks the user to evaluate each fit (excluding references) by eye. If the
% user says the fit is good then the program continues as normal. If the
% user says the fit is bad then an interface to enter new guesses appears
% (the previous guesses are displayed as well). A plot of the new fit with
% the new guess is shown. This loop repeats until the user says the fit is
% good, then the program continues as normal.
% 
% A structure called plotarrays is created for easy of plotting. Each
% array can be specified in the following way:
% plotarrays.fitparameter.(ref or AC)
% 
% The user is asked if they would like to manually input some data to plot
% w0s against. For example power or polarization of AC laser. If the
% incorrect amount of data is entered, the user is given unlimited tries to
% enter the correct amount, feedback is provided. If the program has been
% run before the defualt input for userdata will be the array previously
% entered.
% 
% Plots of (w0s,linewidths,heights,areas) vs. userdata are made with and
% without shaded regions depicting the FWHM of the peaks for the w0s plot.
% 
% The figure and relevant variables are saved.
% ***********************************

close all;
clc;
%% Input parameters
    in_params = inputParser;
    in_params.CaseSensitive = false;
    in_params.addParameter('path','', @ischar);
    in_params.parse(varargin{:});
    
    path = in_params.Results.path;
%% Loading in the data
    % Check the path
    if isempty(path)
        % Ask the user to select the folder of interest '\\ecas.wvu.edu\squol\AC Stark Effect'
        path = uigetdir();
        if path == 0 % User pressed cancel
            cprintf('err', '\nCANCELLED: Folder path selection cancelled.\n');
            return
        end
    end
    
    % Extract and print the name of the folder the user has selected
    split = strsplit(path, '\');
    folder = split(end);
    fprintf(1, ['\nRPLE_Analysis: ', folder{1}, '\n']);
    
    % Load the data
    [data, err] = loadRPLEdata(path);
    if err
        return
    end
    
    % Define total number of scans as N
    N = length(data);
    
    fprintf(1, ['\n' num2str(N) ' RPLE spectra loaded in:']);
    
    % Extract the parameter values for each scan
    data = extractParams(data);
    
%% Determine type of scan
    % scanType: 0 means reference, 1 means AC, 2 means magnet,
    % 3 means AC magnet, and 4 means not AC Stark data
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
    % All other scan types
    Nother = sum(ismember(scanType, 6));
    
    % Print number of scans, and determine type of scans
    % First some assertions to check everything makes sense
    if (Nref + NAC + Nmag + NACmag) ~= 0 && (Nres + NHeNe) ~= 0 && Nother ~= 0
        cprintf('err', '\nERROR: All three types of scans found in the same folder...yikes.\n');
        beep; return
    elseif (Nref + NAC + Nmag + NACmag) ~= 0 && (Nres + NHeNe) ~= 0
        cprintf('err', '\nERROR: AC Stark and excitation laser power scans found in the same folder.\n');
        beep; return
    elseif (Nref + NAC + Nmag + NACmag) ~= 0 && Nother ~= 0
        cprintf('err', '\nERROR: AC Stark and "other" scans found in the same folder.\n');
        beep; return
    elseif (Nres + NHeNe) ~= 0 && Nother ~= 0
        cprintf('err', '\nERROR: Excitation laser power and "other" scans found in the same folder.\n');
        beep; return
    elseif (Nref + NAC + Nmag + NACmag) ~= 0
        ACStarkScans = 1;
        excitationLaserPowerScans = 0;
        otherScans = 0;
        fprintf(1, '\tAC Stark scans\n');
        fprintf(1, ['\t\t\t\t\t\t\t' num2str(Nref) ' reference\n']);
        fprintf(1, ['\t\t\t\t\t\t\t' num2str(NAC) ' AC\n']);
        fprintf(1, ['\t\t\t\t\t\t\t' num2str(Nmag) ' magnet\n']);
        fprintf(1, ['\t\t\t\t\t\t\t' num2str(NACmag) ' AC magnet\n']);
    elseif (Nres + NHeNe) ~= 0
        ACStarkScans = 0;
        excitationLaserPowerScans = 1;
        otherScans = 0;
        fprintf(1, '\texcitation laser power scans\n');
        fprintf(1, ['\t\t\t\t\t\t\t' num2str(Nres) ' resonant OD\n']);
        fprintf(1, ['\t\t\t\t\t\t\t' num2str(NHeNe) ' HeNe OD\n']);
    elseif Nother ~= 0
        ACStarkScans = 0;
        excitationLaserPowerScans = 0;
        otherScans = 1;
        fprintf(1, '\tother scans\n');
        fprintf(1, ['\t\t\t\t\t\t\t' num2str(Nother) ' total\n']);
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
    answer0 = questdlg('Would you like to visually check each fit?',...
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
    end
    
    % For loop to fit all of the spectra
    for i = 1:N
        
        % Data for fit
        xData = data(i).x;
        yData = data(i).y;
        yErr = data(i).sy;
        
        % Fit with a single Lorentzian
        [f, gof, guess] = fitLorentzian(xData, yData, yErr);
        
        % Fit R and L detection with 1 peak only if AC Stark data
        if ACStarkScans && ((data(i).detPol == 'R') || (data(i).detPol == 'L'))
            tolerance = 0.6;
        end
        
        % Evaluate rsquared to see if need to fit with 2 Lorentzians
        badfit = 0;
        if gof.rsquare <= tolerance
            
            % Fit the data with a sum of two Lorentzians
            [f, gof, guess] = fit2Lorentzians(xData, yData, yErr);
            
            % Record if fit is still bad
            if gof.rsquare <= tolerance
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
            [numPeaks(i), cancelled] = manualFitting(xData, yData, yErr, guess, numPeaks(i));
            if cancelled
                return
            end
        end
        
        % Assign parameters to fitvalues structure
        fitvalues(i) = assignFitValues(data(i), fitvalues(i), f, numPeaks(i));
    end
    
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
        %{
        if Nres ~= 0 && NHeNe == 0 % Resonant OD plots
            plots_res(plotarrays)
        elseif NHeNe == 0 && Nres == 0 % HeNe OD plots
            plots_HeNe(plotarrays)
        else % Resonant and HeNe plots
            plots_resHeNe(plotarrays)
        end
        %}
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
    savefig(h, [ path '\' folder{1} ' Plots' '.fig']);
    
    % Save variables from the workspace
    if ACStarkScans
        save([ path '\' folder{1}, ' RPLE data', '.mat'], 'data', 'fitvalues', 'plotarrays',...
            'Nref', 'NAC', 'Nmag', 'NACmag', 'scanType', 'numPeaks', 'tolerance', 'refmean');
    else
        save([ path '\' folder{1}, ' RPLE data', '.mat'], 'data', 'fitvalues', 'plotarrays',...
            'Nother', 'scanType', 'numPeaks', 'tolerance');
    end
    
    fprintf(1, '\nData and figures saved!\n');
    
end