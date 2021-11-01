function [f, gof, guess] = fitLorentzian(xData, yData, yErr, varargin)
%% fitLorentzian.m
% Function to fit RPLE data with a single Lorentzian. If none are supplied,
% the initial guesses are calculated based on the data.

%% Input Parameters
    p = inputParser;
    p.CaseSensitive = false;
    p.addParameter('guess',[], @isrow);
    p.parse(varargin{:});
    
    guess = p.Results.guess;
    
%% 
    % Define the Lorentzian function
    Lorentz = @(A,gamma,x0,B,x)...
        (A/pi)*(0.5*gamma)*(((x-x0).^2+0.25*(gamma^2)).^-1)+B;
    
    if isempty(guess) % No guesses supplied, so calculate
        % Initial guesses
        [~, maxIntIndex] = max(yData); % x index for max
        x0 = xData(maxIntIndex); % Peak frequency
        B = min(yData); % Background

        % Linewidth
        halfmax = 0.5*(max(yData)+min(yData)); % Half height
        index1 = find(yData >= halfmax, 1, 'first'); % First time above half height
        index2 = find(yData >= halfmax, 1, 'last'); % Last time above half height
        gamma = xData(index2)-xData(index1); % Approximate linewidth

        A = 0.5*pi*gamma*(max(yData)-B); % Solve the equation for this last guess
    elseif length(guess) ~= 4
        cprintf('err', ['\nERROR: 4 parameters needed to fit, only ' num2str(length(guess)) ' supplied.\n']);
        beep; return
    else
        A = guess(1);
        gamma = guess(2);
        x0 = guess(3);
        B = guess(4);
    end
    
    % Fit with a single Lorentzian.
    [f, gof] = fit(xData,yData,Lorentz,'StartPoint',[A, gamma, x0, B],'Weights',max(yErr)./yErr,...
        'Lower',[0 0 0 0]);
    
    % Ouput the guesses
    guess = [A, gamma, x0, B];
    
end