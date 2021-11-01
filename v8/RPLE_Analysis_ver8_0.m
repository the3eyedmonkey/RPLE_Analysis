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
% ******* ver7.0 (11/2/21) **********
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
            cprintf('err', '\nERROR: Please select the folder DIRECTLY ABOVE the folders with data in them.\n');
            beep; return
        end
    end
    
    % Extract and print the name of the folder the user has selected
    split = strsplit(path, '\');
    folder = split(end);
    fprintf(1, ['\nRPLE_Analysis: ', folder{1}, '\n']);
    
    % Load the data
    data = loadRPLEdata(path);
    
    % Define total number of scans as N
    N = length(data);
    
    fprintf(1, ['\n' num2str(N) ' RPLE spectra loaded in:\n']);
    
%% Extract the parameter values for each scan
    data = extractParams(data);
    
%% Determine type of scan
    % scanType: 0 means reference, 1 means AC, 2 means magnet, 3 means AC
    % magnet, and 4 means not AC Stark data.
    scanType = determineScanType(data);
    
    % Determine the number of reference and AC scans.
    Nref = sum(ismember(scanType,0));
    NAC = sum(ismember(scanType,1));
    Nmag = sum(ismember(scanType,2));
    NACmag = sum(ismember(scanType,3));
    Nother = sum(ismember(scanType,4));
    
    if Nother ~= 0
        fprintf(1, ['\t\t\t\t\t\t\t', num2str(Nother), ' non AC Stark scans\n']);
    else
        fprintf(1, ['\t\t\t\t\t\t\t', num2str(Nref), ' reference\n']);
        fprintf(1, ['\t\t\t\t\t\t\t', num2str(NAC), ' AC\n']);
        fprintf(1, ['\t\t\t\t\t\t\t', num2str(Nmag), ' magnet\n']);
        fprintf(1, ['\t\t\t\t\t\t\t', num2str(NACmag), ' AC magnet\n']);
    end
    
%% Determine value to use for fitting tolerance
    % If the program has been ran before provide option to choose the same.
    if ranBefore(path)
        [~, usertolerance] = ranBefore(path);
        answertol = questdlg(sprintf([ 'What value for user tolerance would you like to use?'...
            '\n\n0.92 is default'...
            '\n0.6 to fit with one peak (Linear AC or only R/L detection)'...
            '\n' num2str(usertolerance) ' is previously entered value' ]),...
            'User tolerance value entry',...
            '0.92', num2str(usertolerance), 'Custom', num2str(usertolerance));
    else
        answertol = questdlg(sprintf([ 'What value for user tolerance would you like to use?'...
            '\n\n0.92 is default'...
            '\n0.6 to fit with one peak (Linear AC or only R/L detection)' ]),...
            'User tolerance value entry',...
            '0.92', '0.60', 'Custom', '0.92');
    end
    
    % Assign the value the user selected.
    switch answertol
        case '0.92'
            usertolerance = 0.92;
        case '0.60'
            usertolerance = 0.6;
        case '' %User closed the dialog box
            cprintf('err', '\nERROR: Please enter a user tolerance value.\n');
            beep; return
        case 'Custom' %Let the user enter their desired value manually
            answercustomtol = inputdlg('User tolerance value:',...
                'Custom user tolerance value entry',...
                [1 50]);
            usertolerance = str2double(answercustomtol{1});
    end
    
%% Fitting with Lorentzian(s) and extracting parameters
    % Array for feedback on number of peaks in each spectrum
    numPeaks = zeros(1, N);
    
    % Create an array called 'fitvalues' for all of the fit parameters
    chars = char(1:N);
    fitvalues = struct('filename',cellstr((chars(1:N))'),...
            'w0s',cellstr((chars(1:N))'),...
            'linewidths',cellstr((chars(1:N))'),...
            'heights',cellstr((chars(1:N))'),...
            'areas',cellstr((chars(1:N))'),...
            'B',cellstr((chars(1:N))'));
    
    % Define a global feedback variable to know if any scans are fit with
    % the sum of two Lorentzians. And if any scans could not be fit well.
    feedbackglobal = 0;
    
    fprintf(1, '\nFitting progress:\n');
    
    % Ask the user whether this is important data that needs each fit to be
    % checked manually by the user. (Basically only necessary if this is
    % final data for a paper)
    answer0 = questdlg('Would you like to fit this data manually?',...
        'Manual fitting');
    
    switch answer0
        case 'No'
            manualfit = 0;
        case 'Cancel'
            cprintf('err', '\nERROR: Please select either yes or no when asked if you would like to fit scans manually.\n');
            beep; return
        case '' % User closed the dialog box
            cprintf('err', '\nERROR: Please select either yes or no when asked if you would like to fit scans manually.\n');
            beep; return
        case 'Yes' % This is important data the user would like to check each fit manually
            manualfit = 1;
    end
    
    % For loop to fit all of the spectra
    for i = 1:N
        fitvalues(i).filename = data(i).filename;
        
        % Data for fit
        xData = data(i).x;
        yData = data(i).y;
        yErr = data(i).sy;
        
        % Fit with a single Lorentzian
        [f, gof, guess] = fitLorentzian(xData, yData, yErr);
        
        % Control fits based on detection polarization
        if (data(i).detPol == 'X') || (data(i).detPol == 'Y') % Fit X and Y detection with 2 peaks
            tolerance = usertolerance;
        elseif (data(i).detPol == 'R') || (data(i).detPol == 'L') % Fit R and L detection with 1 peak
            tolerance = 0.6;
        else
            tolerance = usertolerance;
        end
        
        % Evaluate rsquared to see if need to fit with 2 Lorentzians
        badfit = 0;
        if any(scanType(i) == [1 2 3]) && gof.rsquare <= tolerance
            feedbackglobal = 1; % Global knowledge of if any scans are fit with 2 Lorentzians
            
            % Fit the data with a sum of two Lorentzians
            [f, gof, guess] = fit2Lorentzians(xData, yData, yErr);
            
            % Check to see if the fit is good
            if gof.rsquare <= tolerance
                badfit = 1;
                fprintf(1, ['\t\t\t\t' data(i).filename ' has not been fit well\n']);
            end
            
            % Set the number of peaks
            numPeaks(i) = 2;
        else
            numPeaks(i) = 1;
        end
        
        % Plot the raw data with fits
        figure('Name', data(i).filename,'WindowStyle', 'docked', 'numbertitle', 'off');
        hold on
        h0 = errorbar(xData,yData,yErr,'b.','Capsize',0.1);
        h1 = plot(f);
        
        % Define the Lorentzian function for plotting
        Lorentz = @(A,gamma,x0,B,x)...
            (A/pi)*(0.5*gamma)*(((x-x0).^2+0.25*(gamma^2)).^-1)+B;
        
        % Define the sum of two Lorentzians for plotting
        Lorentz2 = @(A1,A2,gamma1,gamma2,x01,x02,B,x)...
            (A1/pi)*(0.5*gamma1)*(((x-x01).^2+0.25*(gamma1^2)).^-1)+...
            (A2/pi)*(0.5*gamma2)*(((x-x02).^2+0.25*(gamma2^2)).^-1)+B;
        
        % If the fit is bad, show the guessed Lorentzians
        if badfit == 1 && length(q) == 1
            h2 = plot(xData,...
                Lorentz2(guess(1),guess(2),guess(3),guess(4),guess(5),guess(6),guess(7),xData),...
                'k'); % Guess will appear black
            legend([h0,h1,h2],'data', 'fit', 'guess')
        elseif badfit == 1 && length(q) == 2
            h2 = plot(xData,...
                Lorentz2(guess(1),guess(2),guess(3),guess(4),guess(5),guess(6),guess(7),xData),...
                'g'); % Guess will appear green
            legend([h0,h1,h2],'data', 'fit', 'guess')
        else
            legend([h0,h1],'data','fit')
        end
        
        xlabel('Frequency (GHz)');
        ylabel('Intensity (arb. units)');
        
        % If the user has specified that this is important data, allow them
        % to look at the fit and decide if it is good enough
        if manualfit == 1 && any(scanType(i) == [1 2 3 4]) % Only worry about this for non-reference scans
            answerfit = questdlg('Is this fit okay?','Manual Fitting');
            
            switch answerfit
                case 'No' % The fit is bad and the user would like to enter guesses manually
                    manualbadfit = 1;
                    % Variables to tell how the user would like to refit
                    singletosingle = 0;
                    singletodouble = 0;
                    doubletosingle = 0;
                    doubletodouble = 0;
                    
                    % Ask how many Lorentzians to fit with
                    answerdouble = questdlg('How many Lorentzians would you like to fit with?',...
                        'Manual Fitting',...
                        'One','Two','Two');

                    switch answerdouble
                        case 'One' % This seems unlikely, so I will table it unless it comes up
                            if numPeaks(i) == 1
                                singletosingle = 1;
                            elseif numPeaks(i) == 2
                                doubletosingle = 1;
                                numPeaks(i) = 1;
                            end
                        case ''
                            cprintf('err', '\nERROR: Please select either yes or no when asked how many Lorentzians to fit with.\n');
                            beep; return
                        case 'Two'
                            if numPeaks(i) == 1
                                singletodouble = 1;
                                numPeaks(i) = 2;
                            elseif numPeaks(i) == 2
                                doubletodouble = 1;
                            end
                    end   
                case 'Cancel'
                    cprintf('err', '\nERROR: Please select either yes or no when asked if the fit is okay.\n');
                    beep; return
                case '' % User closed the dialog box
                    cprintf('err', '\nERROR: Please select either yes or no when asked if the fit is okay.\n');
                    beep; return
                case 'Yes' % The fit is good, so move on
                    manualbadfit = 0;
            end
            
            % Now we know how many Lorentzians the initial data was fit
            % with, and how many the user would like to manually fit with.
            % Next run a while loop until to refit until the user is
            % satisfied with the results.
            
            while manualbadfit == 1 % Repeat until the user says the fit is good enough
                answerboxdim = 50;
                if singletosingle == 1 || doubletosingle == 1 %The user wants to fit with a single Lorentzian
                    % Previous guesses
                    if length(guess) == 7 % If previously fit with 2 Lorentzians
                        A = guess(1);
                        gamma = guess(3);
                        x0 = guess(5);
                        guess(6) = [];
                        guess(4) = [];
                        guess(2) = [];
                    else
                        A = guess(1);
                        gamma = guess(2);
                        x0 = guess(3);
                    end
                    
                    % User input
                    prompt = {[ 'Guess for A (Previous was ' num2str(A) '):' ],...
                            [ 'Guess for gamma (Previous was ' num2str(gamma) '):' ],...
                            [ 'Guess for w0 (Previous was ' num2str(x0) '):' ]};
                    title = 'Manual Guess Entry';
                    dims = [1 answerboxdim; 1 answerboxdim; 1 answerboxdim];
                    definput = {num2str(A),num2str(gamma),num2str(x0)};
                    opts.Resize = 'on';
                    opts.WindowStyle = 'normal';
                    
                    answersinglefit = inputdlg(prompt,title,dims,definput,opts);
                    
                    % User guesses
                    for j = 1:length(answersinglefit)
                        guess(j) = str2double(answersinglefit{j});
                    end
                    guess(4) = min(yData);
                    
                    % Refit with user guesses
                    [f, ~, guess] = fitLorentzian(xData, yData, yErr, 'guess', guess);
                    
                    % Plot the fit again so the user can look at it.
                    clf % Clear the current figure window
                    hold on
                    h0 = errorbar(xData,yData,yErr,'b.','Capsize',0.1); % Overwrite old figure
                    h1 = plot(f);
                    
                    % Plot the guess as well to aid the user in fitting
                    h2 = plot(xData, Lorentz(guess(1),guess(2),guess(3),guess(4),xData), 'g');
                    legend([h0,h1,h2],'data', 'fit', 'guess')
                    
                    % Ask again if the fit is good.
                    answerfit = questdlg('Is this fit okay?','Manual fitting');
                    
                    switch answerfit
                        case 'No'
                            manualbadfit = 1;
                        case 'Cancel'
                            cprintf('err', '\nERROR: Please select either yes or no when asked if the fit is okay.\n');
                            beep; return
                        case ''
                            cprintf('err', '\nERROR: Please select either yes or no when asked if the fit is okay.\n');
                            beep; return
                        case 'Yes'
                            numPeaks(i) = 1; % Keep track of number of peaks for analysis later on
                            
                            % Get rid of the guess on the plot
                            clf
                            h0 = errorbar(xData,yData,yErr,'b.','Capsize',0.1);
                            hold on
                            h1 = plot(f);
                            legend([h0,h1],'data', 'fit')
                            
                            manualbadfit = 0; % Exit the while loop
                    end
                elseif singletodouble == 1 || doubletodouble == 1 %The user wants to fit with double Lorentzian
                    % Previous guesses
                    if length(guess) == 4 % If previously fit with single Lorentzian
                        A1 = guess(1);
                        A2 = 0;
                        gamma1 = guess(2);
                        gamma2 = 0;
                        x01 = guess(3);
                        x02 = 0;
                        B = guess(4);
                    else
                        A1 = guess(1);
                        A2 = guess(2);
                        gamma1 = guess(3);
                        gamma2 = guess(4);
                        x01 = guess(5);
                        x02 = guess(6);
                        B = guess(7);
                    end
                    
                    % User input
                    prompt = {[ 'Guess for A1 (Previous was ' num2str(A1) '):' ],...
                            [ 'Guess for A2 (Previous was ' num2str(A2) '):' ],...
                            [ 'Guess for gamma1 (Previous was ' num2str(gamma1) '):' ],...
                            [ 'Guess for gamma2 (Previous was ' num2str(gamma2) '):' ],...
                            [ 'Guess for w01 (Previous was ' num2str(x01) '):' ],...
                            [ 'Guess for w02 (Previous was ' num2str(x02) '):' ],...
                            [ 'Guess for B (Previous was ' num2str(B) '):' ]};
                    title = 'Manual Guess Entry';
                    dims = [1 answerboxdim; 1 answerboxdim; 1 answerboxdim; 1 answerboxdim;...
                        1 answerboxdim; 1 answerboxdim; 1 answerboxdim];
                    definput = {num2str(A1),num2str(A2),num2str(gamma1),num2str(gamma2),...
                            num2str(x01),num2str(x02),num2str(B)};
                    opts.Resize = 'on';
                    opts.WindowStyle = 'normal';
                    
                    answerdoublefit = inputdlg(prompt,title,dims,definput,opts);
                    
                    % User guesses
                    for j = 1:length(answerdoublefit)
                        guess(j) = str2double(answerdoublefit{j});
                    end
                    
                    % Refit with user guesses
                    [f, ~, guess] = fit2Lorentzians(xData, yData, yErr, 'guess', guess);
                    
                    % Plot the fit again so the user can look at it.
                    clf % Clear the current figure window
                    hold on
                    h0 = errorbar(xData,yData,yErr,'b.','Capsize',0.1); % Overwrite old figure
                    h1 = plot(f);
                    
                    % Plot the guess as well to aid the user in fitting
                    h2 = plot(xData,...
                        Lorentz2(guess(1),guess(2),guess(3),guess(4),guess(5),guess(6),guess(7),xData),...
                        'g');
                    legend([h0,h1,h2],'data', 'fit', 'guess')
                    
                    % Ask again if the fit is good.
                    answerfit = questdlg('Is this fit okay?','Manual fitting');
                    
                    switch answerfit
                        case 'No'
                            manualbadfit = 1;
                        case 'Cancel'
                            cprintf('err', '\nERROR: Please select either yes or no when asked if the fit is okay.\n');
                            beep; return
                        case ''
                            cprintf('err', '\nERROR: Please select either yes or no when asked if the fit is okay.\n');
                            beep; return
                        case 'Yes'
                            numPeaks(i) = 1; % Keep track of number of peaks for analysis later on
                            
                            % Get rid of the guess on the plot
                            clf
                            h0 = errorbar(xData,yData,yErr,'b.','Capsize',0.1);
                            hold on
                            h1 = plot(f);
                            legend([h0,h1],'data', 'fit')
                            
                            manualbadfit = 0; % Exit the while loop
                    end
                end
            end
        end
        
        % Now output all of the fit parameters into 'fitvalues' structure.
        % If a single Lorentzian is used saves like [value, confint1, confint2]
        % If a double Lorentzian is used saves like
        % [value1, confint11, confint12, value2, confint21, confint22]
        values = coeffvalues(f);
        intervals = confint(f);
        
        if numPeaks(i) == 1
            % Working out height from fitting parameters.
            h = (2*values(1))/(pi*values(2)) + values(4);
            hint1 = (2*intervals(1,1))/(pi*values(2)) + values(4);
            hint2 = (2*intervals(2,1))/(pi*values(2)) + values(4);
            
            % Place fit parameters in the fitvalues structure.
            fitvalues(i).heights = [ h hint1 hint2 ];
            fitvalues(i).areas = [ values(1) intervals(1,1) intervals(2,1) ];
            fitvalues(i).linewidths = [ values(2) intervals(1,2) intervals(2,2) ];
            fitvalues(i).w0s = [ values(3) intervals(1,3) intervals(2,3) ];
            fitvalues(i).B = [ values(4) intervals(1,4) intervals(2,4) ];
        elseif numPeaks(i) == 2
            % Working out height from fitting parameters.
            h1 = (2*values(1))/(pi*values(3)) + (values(2)*values(4))/(2*pi*((values(5)-values(6))^2+(0.5*values(4))^2)) + values(7);
            h1int1 = (2*intervals(1,1))/(pi*values(3)) + (values(2)*values(4))/(2*pi*((values(5)-values(6))^2+(0.5*values(4))^2)) + values(7);
            h1int2 = (2*intervals(2,1))/(pi*values(3)) + (values(2)*values(4))/(2*pi*((values(5)-values(6))^2+(0.5*values(4))^2)) + values(7);
            
            h2 = (values(1)*values(3))/(2*pi*((values(6)-values(5))^2+(0.5*values(3))^2)) + (2*values(2))/(pi*values(4)) + values(7);
            h2int1 = (values(1)*values(3))/(2*pi*((values(6)-values(5))^2+(0.5*values(3))^2)) + (2*intervals(1,2))/(pi*values(4)) + values(7);
            h2int2 = (values(1)*values(3))/(2*pi*((values(6)-values(5))^2+(0.5*values(3))^2)) + (2*intervals(2,2))/(pi*values(4)) + values(7);
            
            
            fitvalues(i).heights = [ h1 h1int1 h1int2...
                h2 h2int1 h2int2 ];
            
            fitvalues(i).areas = [ values(1) intervals(1,1) intervals(2,1)...
                values(2) intervals(1,2) intervals(2,2) ];
            
            fitvalues(i).linewidths = [ values(3) intervals(1,3) intervals(2,3)...
                values(4) intervals(1,4) intervals(2,4) ];
            
            fitvalues(i).w0s = [ values(5) intervals(1,5) intervals(2,5)...
                values(6) intervals(1,6) intervals(2,6) ];
            
            fitvalues(i).B = [ values(7) intervals(1,7) intervals(2,7) ];
        end
    end
    
    
    % Delete parts of 'fitvalues' that were unused (do I still need this?)
    for i = length(fitvalues):-1:1
        if ischar(fitvalues(i).w0s) == 1
            fitvalues(i)=[];
        end
    end
    
    % Check to confirm that data were loaded and fit appropriately (do I
    % still need this)
    if N ~= length(fitvalues)
        cprintf('err', '\nERROR: The length of the data array does not match the length of the fitvalues array.\n');
        beep; return
    end
    
%% Plot all RPLE Spectra
    figure('Name', 'All Scans', 'WindowStyle', 'docked', 'numbertitle', 'off');
    hold on
    
    cmap = jet(N);
    
    % Plot raw data.
    for i = 1:N
        plot(data(i).x,data(i).y,'Color',cmap(i,:))
    end
    
    legend(data(1:end).filename) % Legend labels spectra with filenames
    
    xlabel('Frequency (GHz)');
    ylabel('Intensity (arb. units)');
    hold off
    
%% Make structure of arrays for fit parameters
    % The fitvalues structure can tell you all things about one scan, now
    % we create a plotarrays structure that will tell you one thing about
    % all scans.
    
    plotarrays = struct('indices',struct('ref',chars(1),'AC',chars(1),'mag',chars(1),'ACmag',chars(1)),...
        'ACvalues',struct('AC',chars(1),'ACmag',chars(1)),...
        'magvalues',struct('mag',chars(1),'ACmag',chars(1)),...
        'w0s',struct('ref',chars(1),'referr',chars(1),'AC',chars(1),'ACerr',chars(1),...
            'mag',chars(1),'magerr',chars(1),'ACmag',chars(1),'ACmagerr',chars(1)),...
        'linewidths',struct('ref',chars(1),'referr',chars(1),'AC',chars(1),'ACerr',chars(1),...
            'mag',chars(1),'magerr',chars(1),'ACmag',chars(1),'ACmagerr',chars(1)),...
        'heights',struct('ref',chars(1),'referr',chars(1),'AC',chars(1),'ACerr',chars(1),...
            'mag',chars(1),'magerr',chars(1),'ACmag',chars(1),'ACmagerr',chars(1)),...
        'areas',struct('ref',chars(1),'referr',chars(1),'AC',chars(1),'ACerr',chars(1),...
            'mag',chars(1),'magerr',chars(1),'ACmag',chars(1),'ACmagerr',chars(1)),...
        'B',struct('ref',chars(1),'referr',chars(1),'AC',chars(1),'ACerr',chars(1),...
            'mag',chars(1),'magerr',chars(1),'ACmag',chars(1),'ACmagerr',chars(1)));
    
    % Initialize arrays.
    plotarrays.indices.ref = zeros(1, N);
    plotarrays.indices.AC = zeros(1, 2*N);
    plotarrays.indices.mag = zeros(1, 2*N);
    plotarrays.indices.ACmag = zeros(1, 2*N);
    
    plotarrays.ACvalues.AC = zeros(1, 2*N);
    plotarrays.ACvalues.ACmag = zeros(1, 2*N);
    
    plotarrays.magvalues.mag = zeros(1, 2*N);
    plotarrays.magvalues.ACmag = zeros(1, 2*N);
    
    
    % Reference
    plotarrays.w0s.ref = zeros(1, N);
    plotarrays.linewidths.ref = zeros(1, N);
    plotarrays.heights.ref = zeros(1, N);
    plotarrays.areas.ref = zeros(1, N);
    plotarrays.B.ref = zeros(1, N);
    
    % Reference errobars
    plotarrays.w0s.referr = zeros(1, N);
    plotarrays.linewidths.referr = zeros(1, N);
    plotarrays.heights.referr = zeros(1, N);
    plotarrays.areas.referr = zeros(1, N);
    plotarrays.B.referr = zeros(1, N);
    
    % AC
    plotarrays.w0s.AC = zeros(1, 2*N);
    plotarrays.linewidths.AC = zeros(1, 2*N);
    plotarrays.heights.AC = zeros(1, 2*N);
    plotarrays.areas.AC = zeros(1, 2*N);
    plotarrays.B.AC = zeros(1, 2*N);
    
    % AC error bars
    plotarrays.w0s.ACerr = zeros(1, 2*N);
    plotarrays.linewidths.ACerr = zeros(1, 2*N);
    plotarrays.heights.ACerr = zeros(1, 2*N);
    plotarrays.areas.ACerr = zeros(1, 2*N);
    plotarrays.B.ACerr = zeros(1, 2*N);
    
    % Magnet
    plotarrays.w0s.mag = zeros(1, 2*N);
    plotarrays.linewidths.mag = zeros(1, 2*N);
    plotarrays.heights.mag = zeros(1, 2*N);
    plotarrays.areas.mag = zeros(1, 2*N);
    plotarrays.B.mag = zeros(1, 2*N);
    
    % Magnet error bars
    plotarrays.w0s.magerr = zeros(1, 2*N);
    plotarrays.linewidths.magerr = zeros(1, 2*N);
    plotarrays.heights.magerr = zeros(1, 2*N);
    plotarrays.areas.magerr = zeros(1, 2*N);
    plotarrays.B.magerr = zeros(1, 2*N);
    
    % ACmagnet
    plotarrays.w0s.ACmag = zeros(1, 2*N);
    plotarrays.linewidths.ACmag = zeros(1, 2*N);
    plotarrays.heights.ACmag = zeros(1, 2*N);
    plotarrays.areas.ACmag = zeros(1, 2*N);
    plotarrays.B.ACmag = zeros(1, 2*N);
    
    % ACmagnet error bars
    plotarrays.w0s.ACmagerr = zeros(1, 2*N);
    plotarrays.linewidths.ACmagerr = zeros(1, 2*N);
    plotarrays.heights.ACmagerr = zeros(1, 2*N);
    plotarrays.areas.ACmagerr = zeros(1, 2*N);
    plotarrays.B.ACmagerr = zeros(1, 2*N);
    
    % For loop to run through all scans
    for i = 1:N
        if scanType(i) == 0 %Reference scans
            plotarrays.indices.ref(i) = i;
            
            plotarrays.w0s.ref(i) = fitvalues(i).w0s(1);
            plotarrays.linewidths.ref(i) = fitvalues(i).linewidths(1);
            plotarrays.heights.ref(i) = fitvalues(i).heights(1);
            plotarrays.areas.ref(i) = fitvalues(i).areas(1);
            plotarrays.B.ref(i) = fitvalues(i).B(1);
            
            plotarrays.w0s.referr(i) = fitvalues(i).w0s(3)-fitvalues(i).w0s(2); %confint are symmetric so upper-lower to get error bar
            plotarrays.linewidths.referr(i) = fitvalues(i).linewidths(3)-fitvalues(i).linewidths(2);
            plotarrays.heights.referr(i) = fitvalues(i).heights(3)-fitvalues(i).heights(2);
            plotarrays.areas.referr(i) = fitvalues(i).areas(3)-fitvalues(i).areas(2);
            plotarrays.B.referr(i) = fitvalues(i).B(3)-fitvalues(i).B(2);
            
        elseif scanType(i) == 1 %AC scan
            plotarrays.indices.AC(i) = i;
            plotarrays.ACvalues.AC(i) = data(i).ACvalue;
            
            plotarrays.w0s.AC(i) = fitvalues(i).w0s(1);
            plotarrays.linewidths.AC(i) = fitvalues(i).linewidths(1);
            plotarrays.heights.AC(i) = fitvalues(i).heights(1);
            plotarrays.areas.AC(i) = fitvalues(i).areas(1);
            plotarrays.B.AC(i) = fitvalues(i).B(1);
            
            plotarrays.w0s.ACerr(i) = fitvalues(i).w0s(3)-fitvalues(i).w0s(2);
            plotarrays.linewidths.ACerr(i) = fitvalues(i).linewidths(3)-fitvalues(i).linewidths(2);
            plotarrays.heights.ACerr(i) = fitvalues(i).heights(3)-fitvalues(i).heights(2);
            plotarrays.areas.ACerr(i) = fitvalues(i).areas(3)-fitvalues(i).areas(2);
            plotarrays.B.ACerr(i) = fitvalues(i).B(3)-fitvalues(i).B(2);
            
            if numPeaks(i) == 2 %Fit with 2 Lorentzians
                plotarrays.indices.AC(i + (N-1)) = i;
                plotarrays.ACvalues.AC(i + (N-1)) = data(i).ACvalue;
                
                plotarrays.w0s.AC(i + (N-1)) = fitvalues(i).w0s(4);
                plotarrays.linewidths.AC(i + (N-1)) = fitvalues(i).linewidths(4);
                plotarrays.heights.AC(i + (N-1)) = fitvalues(i).heights(4);
                plotarrays.areas.AC(i + (N-1)) = fitvalues(i).areas(4);
                plotarrays.B.AC(i + (N-1)) = plotarrays.B.AC(i);
                
                plotarrays.w0s.ACerr(i + (N-1)) = fitvalues(i).w0s(6)-fitvalues(i).w0s(5);
                plotarrays.linewidths.ACerr(i + (N-1)) = fitvalues(i).linewidths(6)-fitvalues(i).linewidths(5);
                plotarrays.heights.ACerr(i + (N-1)) = fitvalues(i).heights(6)-fitvalues(i).heights(5);
                plotarrays.areas.ACerr(i + (N-1)) = fitvalues(i).areas(6)-fitvalues(i).areas(5);
                plotarrays.B.ACerr(i + (N-1)) = plotarrays.B.ACerr(i);
            end
            
        elseif scanType(i) == 2 %Magnet scan
            plotarrays.indices.mag(i) = i;
            plotarrays.magvalues.mag(i) = data(i).magvalue;
            
            plotarrays.w0s.mag(i) = fitvalues(i).w0s(1);
            plotarrays.linewidths.mag(i) = fitvalues(i).linewidths(1);
            plotarrays.heights.mag(i) = fitvalues(i).heights(1);
            plotarrays.areas.mag(i) = fitvalues(i).areas(1);
            plotarrays.B.mag(i) = fitvalues(i).B(1);
            
            plotarrays.w0s.magerr(i) = fitvalues(i).w0s(3)-fitvalues(i).w0s(2);
            plotarrays.linewidths.magerr(i) = fitvalues(i).linewidths(3)-fitvalues(i).linewidths(2);
            plotarrays.heights.magerr(i) = fitvalues(i).heights(3)-fitvalues(i).heights(2);
            plotarrays.areas.magerr(i) = fitvalues(i).areas(3)-fitvalues(i).areas(2);
            plotarrays.B.magerr(i) = fitvalues(i).B(3)-fitvalues(i).B(2);
            
            if numPeaks(i) == 2 %Fit with 2 Lorentzians
                plotarrays.indices.mag(i + (N-1)) = i;
                plotarrays.magvalues.mag(i + (N-1)) = data(i).magvalue;
                
                plotarrays.w0s.mag(i + (N-1)) = fitvalues(i).w0s(4);
                plotarrays.linewidths.mag(i + (N-1)) = fitvalues(i).linewidths(4);
                plotarrays.heights.mag(i + (N-1)) = fitvalues(i).heights(4);
                plotarrays.areas.mag(i + (N-1)) = fitvalues(i).areas(4);
                plotarrays.B.mag(i + (N-1)) = plotarrays.B.mag(i);
                
                plotarrays.w0s.magerr(i + (N-1)) = fitvalues(i).w0s(6)-fitvalues(i).w0s(5);
                plotarrays.linewidths.magerr(i + (N-1)) = fitvalues(i).linewidths(6)-fitvalues(i).linewidths(5);
                plotarrays.heights.magerr(i + (N-1)) = fitvalues(i).heights(6)-fitvalues(i).heights(5);
                plotarrays.areas.magerr(i + (N-1)) = fitvalues(i).areas(6)-fitvalues(i).areas(5);
                plotarrays.B.magerr(i + (N-1)) = plotarrays.B.magerr(i);
            end
            
        elseif scanType(i) == 3 %ACmagnet scan
            plotarrays.indices.ACmag(i) = i;
            plotarrays.ACvalues.ACmag(i) = data(i).ACvalue;
            plotarrays.magvalues.ACmag(i) = data(i).magvalue;
            
            plotarrays.w0s.ACmag(i) = fitvalues(i).w0s(1);
            plotarrays.linewidths.ACmag(i) = fitvalues(i).linewidths(1);
            plotarrays.heights.ACmag(i) = fitvalues(i).heights(1);
            plotarrays.areas.ACmag(i) = fitvalues(i).areas(1);
            plotarrays.B.ACmag(i) = fitvalues(i).B(1);
            
            plotarrays.w0s.ACmagerr(i) = fitvalues(i).w0s(3)-fitvalues(i).w0s(2);
            plotarrays.linewidths.ACmagerr(i) = fitvalues(i).linewidths(3)-fitvalues(i).linewidths(2);
            plotarrays.heights.ACmagerr(i) = fitvalues(i).heights(3)-fitvalues(i).heights(2);
            plotarrays.areas.ACmagerr(i) = fitvalues(i).areas(3)-fitvalues(i).areas(2);
            plotarrays.B.ACmagerr(i) = fitvalues(i).B(3)-fitvalues(i).B(2);
            
            if numPeaks(i) == 2 %Fit with 2 Lorentzians
                plotarrays.indices.ACmag(i + (N-1)) = i;
                plotarrays.ACvalues.ACmag(i + (N-1)) = data(i).ACvalue;
                plotarrays.magvalues.ACmag(i + (N-1)) = data(i).magvalue;
                
                plotarrays.w0s.ACmag(i + (N-1)) = fitvalues(i).w0s(4);
                plotarrays.linewidths.ACmag(i + (N-1)) = fitvalues(i).linewidths(4);
                plotarrays.heights.ACmag(i + (N-1)) = fitvalues(i).heights(4);
                plotarrays.areas.ACmag(i + (N-1)) = fitvalues(i).areas(4);
                plotarrays.B.ACmag(i + (N-1)) = plotarrays.B.ACmag(i);
                
                plotarrays.w0s.ACmagerr(i + (N-1)) = fitvalues(i).w0s(6)-fitvalues(i).w0s(5);
                plotarrays.linewidths.ACmagerr(i + (N-1)) = fitvalues(i).linewidths(6)-fitvalues(i).linewidths(5);
                plotarrays.heights.ACmagerr(i + (N-1)) = fitvalues(i).heights(6)-fitvalues(i).heights(5);
                plotarrays.areas.ACmagerr(i + (N-1)) = fitvalues(i).areas(6)-fitvalues(i).areas(5);
                plotarrays.B.ACmagerr(i + (N-1)) = plotarrays.B.ACmagerr(i);
            end            
        end
    end
    
    % Get rid of values that were unassigned.
    for i = N:-1:1 %Reference arrays
        if plotarrays.indices.ref(i) == 0
            plotarrays.indices.ref(i)=[];
        end
        if plotarrays.w0s.ref(i) == 0
            plotarrays.w0s.ref(i)=[];
        end
        if plotarrays.linewidths.ref(i) == 0
            plotarrays.linewidths.ref(i)=[];
        end
        if plotarrays.heights.ref(i) == 0
            plotarrays.heights.ref(i)=[]; 
        end
        if plotarrays.areas.ref(i) == 0
            plotarrays.areas.ref(i)=[]; 
        end
        if plotarrays.B.ref(i) == 0
            plotarrays.B.ref(i)=[]; 
        end
        if plotarrays.w0s.referr(i) == 0
            plotarrays.w0s.referr(i)=[];
        end
        if plotarrays.linewidths.referr(i) == 0
            plotarrays.linewidths.referr(i)=[];
        end
        if plotarrays.heights.referr(i) == 0
            plotarrays.heights.referr(i)=[];
        end
        if plotarrays.areas.referr(i) == 0
            plotarrays.areas.referr(i)=[];
        end
        if plotarrays.B.referr(i) == 0
            plotarrays.B.referr(i)=[];
        end
    end
    
    for i = 2*N:-1:1 %Rest of the arrays
        % AC arrays
        if plotarrays.indices.AC(i) == 0
            plotarrays.indices.AC(i)=[];
        end
        if plotarrays.ACvalues.AC(i) == 0
            plotarrays.ACvalues.AC(i)=[];
        end
        if plotarrays.w0s.AC(i) == 0
            plotarrays.w0s.AC(i)=[];
        end
        if plotarrays.linewidths.AC(i) == 0
            plotarrays.linewidths.AC(i)=[];
        end
        if plotarrays.heights.AC(i) == 0
            plotarrays.heights.AC(i)=[];
        end
        if plotarrays.areas.AC(i) == 0
            plotarrays.areas.AC(i)=[];
        end
        if plotarrays.B.AC(i) == 0
            plotarrays.B.AC(i)=[];
        end
        if plotarrays.w0s.ACerr(i) == 0
            plotarrays.w0s.ACerr(i)=[];
        end
        if plotarrays.linewidths.ACerr(i) == 0
            plotarrays.linewidths.ACerr(i)=[];
        end
        if plotarrays.heights.ACerr(i) == 0
            plotarrays.heights.ACerr(i)=[];
        end
        if plotarrays.areas.ACerr(i) == 0
            plotarrays.areas.ACerr(i)=[];
        end
        if plotarrays.B.ACerr(i) == 0
            plotarrays.B.ACerr(i)=[];
        end
        
        % Magnet arrays
        if plotarrays.indices.mag(i) == 0
            plotarrays.indices.mag(i)=[];
        end
        if plotarrays.magvalues.mag(i) == 0
            plotarrays.magvalues.mag(i)=[];
        end
        if plotarrays.w0s.mag(i) == 0
            plotarrays.w0s.mag(i)=[];
        end
        if plotarrays.linewidths.mag(i) == 0
            plotarrays.linewidths.mag(i)=[];
        end
        if plotarrays.heights.mag(i) == 0
            plotarrays.heights.mag(i)=[];
        end
        if plotarrays.areas.mag(i) == 0
            plotarrays.areas.mag(i)=[];
        end
        if plotarrays.B.mag(i) == 0
            plotarrays.B.mag(i)=[];
        end
        if plotarrays.w0s.magerr(i) == 0
            plotarrays.w0s.magerr(i)=[];
        end
        if plotarrays.linewidths.magerr(i) == 0
            plotarrays.linewidths.magerr(i)=[];
        end
        if plotarrays.heights.magerr(i) == 0
            plotarrays.heights.magerr(i)=[];
        end
        if plotarrays.areas.magerr(i) == 0
            plotarrays.areas.magerr(i)=[];
        end
        if plotarrays.B.magerr(i) == 0
            plotarrays.B.magerr(i)=[];
        end
        
        % ACmagnet arrays
        if plotarrays.indices.ACmag(i) == 0
            plotarrays.indices.ACmag(i)=[];
        end
        if plotarrays.ACvalues.ACmag(i) == 0
            plotarrays.ACvalues.ACmag(i)=[];
        end
        if plotarrays.magvalues.ACmag(i) == 0
            plotarrays.magvalues.ACmag(i)=[];
        end
        if plotarrays.w0s.ACmag(i) == 0
            plotarrays.w0s.ACmag(i)=[];
        end
        if plotarrays.linewidths.ACmag(i) == 0
            plotarrays.linewidths.ACmag(i)=[];
        end
        if plotarrays.heights.ACmag(i) == 0
            plotarrays.heights.ACmag(i)=[];
        end
        if plotarrays.areas.ACmag(i) == 0
            plotarrays.areas.ACmag(i)=[];
        end
        if plotarrays.B.ACmag(i) == 0
            plotarrays.B.ACmag(i)=[];
        end
        if plotarrays.w0s.ACmagerr(i) == 0
            plotarrays.w0s.ACmagerr(i)=[];
        end
        if plotarrays.linewidths.ACmagerr(i) == 0
            plotarrays.linewidths.ACmagerr(i)=[];
        end
        if plotarrays.heights.ACmagerr(i) == 0
            plotarrays.heights.ACmagerr(i)=[];
        end
        if plotarrays.areas.ACmagerr(i) == 0
            plotarrays.areas.ACmagerr(i)=[];
        end
        if plotarrays.B.ACmagerr(i) == 0
            plotarrays.B.ACmagerr(i)=[];
        end
    end
    
    % Check that the reference arrays are correct length.
    if (length(plotarrays.indices.ref) ~= Nref) || (length(plotarrays.indices.ref) ~= length(plotarrays.w0s.ref))
        cprintf('err', '\nERROR: The length of reference plot arrays do not match.\n');
        beep; return
    end

    % Check that the AC arrays are correct length.
    if (length(plotarrays.indices.AC) ~= length(plotarrays.w0s.AC))
        cprintf('err', '\nERROR: The length of AC plot arrays do not match.\n');
        beep; return
    end
    
    % Check that the magnet arrays are the correct length.
    if (length(plotarrays.indices.mag) ~= length(plotarrays.w0s.mag))
        cprintf('err', '\nERROR: The length of magnet plot arrays do not match.\n');
        beep; return
    end
    
    % Check that the ACmagnet arrays are the correct length.
    if (length(plotarrays.indices.ACmag) ~= length(plotarrays.w0s.ACmag))
        cprintf('err', '\nERROR: The length of ACmagnet plot arrays do not match.\n');
        beep; return
    end
    
%% Plot w0s vs. Scan Index
    
    % Deal with the reference scans.
    if Nref == (NAC+Nmag+NACmag) %This only occurs for AOM scans
        smallpwr = 1;
        % Reference to each scan directly before (if magnitude of shift is small).
        refmean = zeros(1, 2*N);
        
        for i = 1:N %iterate through all non reference scans
            if any(scanType(i) == [1 2 3]) && numPeaks(i) == 1
                refmean(i) = fitvalues(i-1).w0s(1);
            elseif any(scanType(i) == [1 2 3]) && numPeaks(i) == 2
                % First peak of nonref scans
                refmean(i) = fitvalues(i-1).w0s(1);
                % Second peak of nonref scans
                refmean(i + (N-1)) = fitvalues(i-1).w0s(1);
            end
        end
        
        % Clean array of values left unassigned.
        for i = 2*N:-1:1
            if refmean(i) == 0
                refmean(i)=[];
            end
        end
    else
        smallpwr = 0;
        % Proceed as usual with the average of all reference scans.
        refmean = mean(plotarrays.w0s.ref);
    end
    
    % Print the refmean to command line.
    fprintf(1, ['\nThe mean of reference scans is ' num2str(refmean) ' GHz .']);
    
    % Make a plot of w0s, with error bars for reference.
    figure('Name', 'w0s vs. Scan Index', 'WindowStyle', 'docked', 'numbertitle', 'off');
    hold on
    
    errorbar(plotarrays.indices.ref,(plotarrays.w0s.ref-mean(refmean)),plotarrays.w0s.referr,'b.','markersize',15,'Capsize',0.1)
    errorbar(plotarrays.indices.AC,(plotarrays.w0s.AC-refmean),plotarrays.w0s.ACerr,'r.','markersize',15,'Capsize',0.1)
    errorbar(plotarrays.indices.mag,(plotarrays.w0s.mag-refmean),plotarrays.w0s.magerr,'g.','markersize',15,'Capsize',0.1)
    errorbar(plotarrays.indices.ACmag,(plotarrays.w0s.ACmag-refmean),plotarrays.w0s.ACmagerr,'m.','markersize',15,'Capsize',0.1)
    hline(0,'k:')
    
    legend('Reference', 'AC', 'magnet', 'ACmagnet');
    
    xlabel('Scan Index');
    ylabel('w0s (GHz)');
    hold off
    
%% Plot (w0s,linewidths,heights,areas) vs. (ACvalues,magvalues)
    Nplots = 0;
    
    % AC scans
    if NAC ~= 0
        Nplots = Nplots + 1;
        % w0s vs. ACvalues
        figure('Name', 'w0s(AC only) vs. AC Stark Laser Power', 'WindowStyle', 'docked', 'numbertitle', 'off');
        hold on

        h1 = errorbar(plotarrays.ACvalues.AC,(plotarrays.w0s.AC-refmean),plotarrays.w0s.ACerr,'r.','markersize',15,'Capsize',0.1);
        h0 = hline(0,'k:');

        legend([h0 h1], 'Reference', 'AC');

        xlabel('AC Stark Laser Power (mW)');
        ylabel('Center Frequency (GHz)');
        hold off
        
        % Linewidths as shaded region (doesn't work)
        %{
        % Make another plot that includes the FWHM(gamma) values of each
        % peak plotted as a shaded area.
        figure('Name', ['w0s vs. AC Stark Laser Power with FWHM'], 'WindowStyle', 'docked', 'numbertitle', 'off');
        
        % Create arrays for +\- gamma/2. These are then plotted as lines,
        % and the area between them is shaded.
        if smallpwr == 0
            upperref = (plotarrays.w0s.ref-refmean)+plotarrays.linewidths.ref/2;
            lowerref = (plotarrays.w0s.ref-refmean)-plotarrays.linewidths.ref/2;
            upperAC = (plotarrays.w0s.AC-refmean)+plotarrays.linewidths.AC/2;
            lowerAC = (plotarrays.w0s.AC-refmean)-plotarrays.linewidths.AC/2;
        elseif smallpwr == 1
            upperAC = (plotarrays.w0s.AC-refmean)+plotarrays.linewidths.AC/2;
            lowerAC = (plotarrays.w0s.AC-refmean)-plotarrays.linewidths.AC/2;
        end
        
        % Calculate the number of AC scans fit with a single Lorentzian.
        nsinglepks = 2*NAC - length(plotarrays.w0s.AC);
        
        % Create the shaded regions for linewidths.
        hold on
        if smallpwr == 0
            h0 = fill([linspace(0,10,Nref) fliplr(linspace(0,10,Nref))], [upperref fliplr(lowerref)], [0.93 0.93 0.93], 'linestyle', 'none');
            fill([plotarrays.w0s.AC(1:NAC)-refmean fliplr(plotarrays.w0s.AC(1:NAC)-refmean)],...
                [upperAC(1:NAC) fliplr(lowerAC(1:NAC))], [0.93 0.93 0.93], 'linestyle', 'none');
            fill([plotarrays.w0s.AC((end/2)+(nsinglepks/2)+1:end)-refmean fliplr(plotarrays.w0s.AC((end/2)+(nsinglepks/2)+1:end))-refmean],...
                [upperAC((end/2)+(nsinglepks/2)+1:end) fliplr(lowerAC((end/2)+(nsinglepks/2)+1:end))], [0.93 0.93 0.93], 'linestyle', 'none')
        elseif smallpwr == 1
            h0 = fill([plotarrays.w0s.AC(1:NAC) fliplr(plotarrays.w0s.AC(1:NAC))],...
                [upperAC(1:NAC) fliplr(lowerAC(1:NAC))], [0.93 0.93 0.93], 'linestyle', 'none');
            fill([plotarrays.w0s.AC((end/2)+(nsinglepks/2)+1:end) fliplr(plotarrays.w0s.AC((end/2)+(nsinglepks/2)+1:end))],...
                [upperAC((end/2)+(nsinglepks/2)+1:end) fliplr(lowerAC((end/2)+(nsinglepks/2)+1:end))], [0.93 0.93 0.93], 'linestyle', 'none')
        end
        
        % Add the w0s to the plot.
        hold all
        h1 = hline(0,'k:');
        h2 = errorbar(plotarrays.ACvalues.AC,(plotarrays.w0s.AC-refmean),plotarrays.w0s.ACerr,'r.','markersize',15,'Capsize',0.1);
        legend([h0 h1 h2], '\gamma  (FWHM)', 'Reference', 'AC');
        
        xlabel('AC Stark Laser Power (mW)');
        ylabel('w0s (GHz)');
        hold off
        %}
        
        % linewidths vs. ACvalues
        figure('Name', 'Linewidths vs. AC Stark Laser Power', 'WindowStyle', 'docked', 'numbertitle', 'off');
        hold on
        h1 = errorbar(plotarrays.ACvalues.AC,plotarrays.linewidths.AC,plotarrays.linewidths.ACerr,'r.','markersize',15,'Capsize',0.1);
        h0 = hline(mean(plotarrays.linewidths.ref),'k:');
        
        legend([h0 h1], 'Reference', 'AC');
        
        xlabel('AC Stark Laser Power (mW)');
        ylabel('Linewidth (GHz)');
        hold off
        
        % heights-background vs. ACvalues
        figure('Name', 'Heights vs. AC Stark Laser Power', 'WindowStyle', 'docked', 'numbertitle', 'off');
        hold on
        h1 = errorbar(plotarrays.ACvalues.AC,plotarrays.heights.AC-plotarrays.B.AC,plotarrays.heights.ACerr,'r.','markersize',15,'Capsize',0.1);
        h0 = hline(mean(plotarrays.heights.ref-plotarrays.B.ref),'k:');
        
        legend([h0 h1], 'Reference', 'AC');
        
        xlabel('AC Stark Laser Power (mW)');
        ylabel('Peak Height (arb. units)');
        hold off
        
        % areas vs. ACvalues
        figure('Name', 'Areas vs. AC Stark Laser Power', 'WindowStyle', 'docked', 'numbertitle', 'off');
        hold on
        h1 = errorbar(plotarrays.ACvalues.AC,plotarrays.areas.AC,plotarrays.areas.ACerr,'r.','markersize',15,'Capsize',0.1);
        h0 = hline(mean(plotarrays.areas.ref),'k:');
        
        legend([h0 h1], 'Reference', 'AC');
        
        xlabel('AC Stark Laser Power (mW)');
        ylabel('Peak Area (Int*GHz)');
        hold off
    end
    
    % Magnet scans
    if Nmag ~= 0
        Nplots = Nplots + 1;
        % w0s vs. magvalues
        figure('Name', 'w0s(magnet only) vs. Magnetic Field Strength', 'WindowStyle', 'docked', 'numbertitle', 'off');
        hold on
        
        h1 = errorbar(plotarrays.magvalues.mag,(plotarrays.w0s.mag-refmean),plotarrays.w0s.magerr,'g.','markersize',15,'Capsize',0.1);
        h0 = hline(0,'k:');
        
        legend([h0 h1], 'Reference', 'magnet');
        
        xlabel('Magnetic Field Strength (mT)');
        ylabel('Center Frequency (GHz)');
        hold off
        
        % linewidths vs. magvalues
        figure('Name', 'Linewidths vs. Magnetic Field Strength', 'WindowStyle', 'docked', 'numbertitle', 'off');
        hold on
        h1 = errorbar(plotarrays.magvalues.mag,plotarrays.linewidths.mag,plotarrays.linewidths.magerr,'g.','markersize',15,'Capsize',0.1);
        h0 = hline(mean(plotarrays.linewidths.ref),'k:');
        
        legend([h0 h1], 'Reference', 'magnet');
        
        xlabel('Magnetic Field Strength (mT)');
        ylabel('Linewidth (GHz)');
        hold off
        
        % heights-background vs. magvalues
        figure('Name', 'Heights vs. Magnetic Field Strength', 'WindowStyle', 'docked', 'numbertitle', 'off');
        hold on
        h1 = errorbar(plotarrays.magvalues.mag,plotarrays.heights.mag-plotarrays.B.mag,plotarrays.heights.magerr,'g.','markersize',15,'Capsize',0.1);
        h0 = hline(mean(plotarrays.heights.ref-plotarrays.B.ref),'k:');
        
        legend([h0 h1], 'Reference', 'magnet');
        
        xlabel('Magnetic Field Strength (mT)');
        ylabel('Peak Height (GHz)');
        hold off
        
        % areas vs. magvalues
        figure('Name', 'Areas vs. Magnetic Field Strength', 'WindowStyle', 'docked', 'numbertitle', 'off');
        hold on
        h1 = errorbar(plotarrays.magvalues.mag,plotarrays.areas.mag,plotarrays.areas.magerr,'g.','markersize',15,'Capsize',0.1);
        h0 = hline(mean(plotarrays.areas.ref),'k:');
        
        legend([h0 h1], 'Reference', 'magnet');
        
        xlabel('Magnetic Field Strength (mT)');
        ylabel('Peak Area (Int*GHz)');
        hold off
    end
    
    % ACmagnet scans
    if NACmag ~= 0
        Nplots = Nplots + 1;
        % w0s
        figure('Name', 'Center Frequencies', 'WindowStyle', 'docked', 'numbertitle', 'off');
        
        % AC scans
        h1 = plot3(plotarrays.ACvalues.AC,zeros(1,length(plotarrays.ACvalues.AC)),...
            (plotarrays.w0s.AC-refmean),'r.','markersize',15);
        hold on
        % Magnet scans
        h2 = plot3(zeros(1,length(plotarrays.magvalues.mag)),plotarrays.magvalues.mag,...
            (plotarrays.w0s.mag-refmean),'g.','markersize',15);
        % ACmagnet scans
        h3 = plot3(plotarrays.ACvalues.ACmag,plotarrays.magvalues.ACmag,...
            (plotarrays.w0s.ACmag-refmean),'m.','markersize',15);

        legend([h1 h2 h3], 'AC', 'magnet', 'ACmagnet');

        xlabel('AC Stark Laser Power (mW)');
        ylabel('Magnetic Field Strength (mT)');
        zlabel('Center Frequency (GHz)');
        grid on
        hold off
        
        
        % linewidths
        figure('Name', 'Linewidths', 'WindowStyle', 'docked', 'numbertitle', 'off');
        
        % AC scans
        h1 = plot3(plotarrays.ACvalues.AC,zeros(1,length(plotarrays.ACvalues.AC)),...
            (plotarrays.linewidths.AC),'r.','markersize',15);
        hold on
        % Magnet scans
        h2 = plot3(zeros(1,length(plotarrays.magvalues.mag)),plotarrays.magvalues.mag,...
            (plotarrays.linewidths.mag),'g.','markersize',15);
        % ACmagnet scans
        h3 = plot3(plotarrays.ACvalues.ACmag,plotarrays.magvalues.ACmag,...
            (plotarrays.linewidths.ACmag),'m.','markersize',15);

        legend([h1 h2 h3], 'AC', 'magnet', 'ACmagnet');

        xlabel('AC Stark Laser Power (mW)');
        ylabel('Magnetic Field Strength (mT)');
        zlabel('Linewidth (GHz)');
        grid on
        hold off
        
        
        % heights-background
        figure('Name', 'Peak Heights', 'WindowStyle', 'docked', 'numbertitle', 'off');
        
        % AC scans
        h1 = plot3(plotarrays.ACvalues.AC,zeros(1,length(plotarrays.ACvalues.AC)),...
            (plotarrays.heights.AC-plotarrays.B.AC),'r.','markersize',15);
        hold on
        % Magnet scans
        h2 = plot3(zeros(1,length(plotarrays.magvalues.mag)),plotarrays.magvalues.mag,...
            (plotarrays.linewidths.mag-plotarrays.B.mag),'g.','markersize',15);
        % ACmagnet scans
        h3 = plot3(plotarrays.ACvalues.ACmag,plotarrays.magvalues.ACmag,...
            (plotarrays.linewidths.ACmag-plotarrays.B.ACmag),'m.','markersize',15);

        legend([h1 h2 h3], 'AC', 'magnet', 'ACmagnet');

        xlabel('AC Stark Laser Power (mW)');
        ylabel('Magnetic Field Strength (mT)');
        zlabel('Peak Height (arb. units)');
        grid on
        hold off
        
        
        % areas
        figure('Name', 'Peak Areas', 'WindowStyle', 'docked', 'numbertitle', 'off');
        
        % AC scans
        h1 = plot3(plotarrays.ACvalues.AC,zeros(1,length(plotarrays.ACvalues.AC)),...
            (plotarrays.areas.AC),'r.','markersize',15);
        hold on
        % Magnet scans
        h2 = plot3(zeros(1,length(plotarrays.magvalues.mag)),plotarrays.magvalues.mag,...
            (plotarrays.areas.mag),'g.','markersize',15);
        % ACmagnet scans
        h3 = plot3(plotarrays.ACvalues.ACmag,plotarrays.magvalues.ACmag,...
            (plotarrays.areas.ACmag),'m.','markersize',15);

        legend([h1 h2 h3], 'AC', 'magnet', 'ACmagnet');

        xlabel('AC Stark Laser Power (mW)');
        ylabel('Magnetic Field Strength (mT)');
        zlabel('Peak Area (Int*GHz)');
        grid on
        hold off
    end

%% Save the figure window and the workspace
    
    % Array with the figure handles.
    Nplots = N + 2 + 4*Nplots;
    h = 1:Nplots;
    
    %Save the figure.
    savefig(h, [ path '\' folder{1} ' Plots' '.fig']);
    
    %Save specific variables from the workspace.
    if feedbackglobal == 0
        save([ path '\' folder{1}, ' RPLE data', '.mat'], 'data', 'fitvalues', 'refmean',...
            'NAC', 'Nref', 'Nmag', 'NACmag', 'scanType', 'numpeaks', 'usertolerance', 'plotarrays');
    elseif feedbackglobal == 1
        save([ path '\' folder{1}, ' RPLE data', '.mat'], 'data', 'fitvalues', 'refmean',...
            'NAC', 'Nref', 'Nmag', 'NACmag', 'scanType', 'numpeaks', 'usertolerance', 'plotarrays');
    end
    
    fprintf(1, '\nData and figures saved!\n');
    
end