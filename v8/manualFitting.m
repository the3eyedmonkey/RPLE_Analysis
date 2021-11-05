function [numpeaks, ft, cancelled] = manualFitting(xData, yData, yErr, guess, numpeaks)
%% manualFitting.m
% Have the user manually supply guesses to fit the data
% Here
% numpeaks = numPeaks(i) from the main function

%% 
    % Let the calling program know if the user manually fit the scan
    ft = 0;
    % Let the calling program know if the user cancelled this
    cancelled = 0;
    
    % Define the Lorentzian function for plotting
    Lorentz = @(A,gamma,x0,B,x)...
        (A/pi)*(0.5*gamma)*(((x-x0).^2+0.25*(gamma^2)).^-1)+B;

    % Define the sum of two Lorentzians for plotting
    Lorentz2 = @(A1,A2,gamma1,gamma2,x01,x02,B,x)...
        (A1/pi)*(0.5*gamma1)*(((x-x01).^2+0.25*(gamma1^2)).^-1)+...
        (A2/pi)*(0.5*gamma2)*(((x-x02).^2+0.25*(gamma2^2)).^-1)+B;
    
    % Ask if the fit is okay
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
                case 'One'
                    if numpeaks == 1
                        singletosingle = 1;
                    elseif numpeaks == 2
                        doubletosingle = 1;
                        numpeaks = 1;
                    end
                case '' % User closed the dialog box
                    cprintf('err', '\nCANCELLED: Manual fitting cancelled.\n');
                    cancelled = 1;
                    return
                case 'Two'
                    if numpeaks == 1
                        singletodouble = 1;
                        numpeaks = 2;
                    elseif numpeaks == 2
                        doubletodouble = 1;
                    end
            end   
        case 'Cancel'
            cprintf('err', '\nCANCELLED: Manual fitting cancelled.\n');
            cancelled = 1;
            return
        case '' % User closed the dialog box
            cprintf('err', '\nCANCELLED: Manual fitting cancelled.\n');
            cancelled = 1;
            return
        case 'Yes' % The fit is good, so move on
            manualbadfit = 0;
    end
    
    % Now we know how many Lorentzians the initial data was fit
    % with, and how many the user would like to manually fit with.
    % Next run a while loop to refit until the user is satisfied.
    
    while manualbadfit == 1 % Repeat until the user says the fit is good enough
        answerboxdim = 50;
        if singletosingle == 1 || doubletosingle == 1 % The user wants to fit with a single Lorentzian
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
                    cprintf('err', '\nCANCELLED: Manual fitting cancelled.\n');
                    cancelled = 1;
                    return
                case ''
                    cprintf('err', '\nCANCELLED: Manual fitting cancelled.\n');
                    cancelled = 1;
                    return
                case 'Yes'
                    numpeaks = 1; % Keep track of number of peaks for analysis later on
                    
                    % Get rid of the guess on the plot
                    clf
                    h0 = errorbar(xData,yData,yErr,'b.','Capsize',0.1);
                    hold on
                    h1 = plot(f);
                    legend([h0,h1],'data', 'fit')
                    
                    ft = 1; % This data fit manually
                    manualbadfit = 0; % Exit the while loop
            end
        elseif singletodouble == 1 || doubletodouble == 1 % The user wants to fit with double Lorentzian
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
                    cprintf('err', '\nCANCELLED: Manual fitting cancelled.\n');
                    cancelled = 1;
                    return
                case ''
                    cprintf('err', '\nCANCELLED: Manual fitting cancelled.\n');
                    cancelled = 1;
                    return
                case 'Yes'
                    numpeaks = 2; % Keep track of number of peaks for analysis later on
                    
                    % Get rid of the guess on the plot
                    clf
                    h0 = errorbar(xData,yData,yErr,'b.','Capsize',0.1);
                    hold on
                    h1 = plot(f);
                    legend([h0,h1],'data', 'fit')
                    
                    ft = 1; % This data fit manually
                    manualbadfit = 0; % Exit the while loop
            end
        end
    end
    
end