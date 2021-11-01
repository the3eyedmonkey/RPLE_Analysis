function [f, gof, guess] = fit2Lorentzians(xData, yData, yErr, varargin)
%% fit2Lorentzians.m
% Function to fit RPLE data with a sum of two Lorentzians. If none are
% supplied, the initial guesses are calculated based on the data.

%% Input Parameters
    p = inputParser;
    p.CaseSensitive = false;
    p.addParameter('guess',[], @isrow);
    p.parse(varargin{:});
    
    guess = p.Results.guess;
    
%% 
    % Define the sum of two Lorentzians
    Lorentz2 = @(A1,A2,gamma1,gamma2,x01,x02,B,x)...
        (A1/pi)*(0.5*gamma1)*(((x-x01).^2+0.25*(gamma1^2)).^-1)+...
        (A2/pi)*(0.5*gamma2)*(((x-x02).^2+0.25*(gamma2^2)).^-1)+B;
    
    if isempty(guess) % No guesses supplied, so calculate
        % Use the findpeaks() function to search for peaks in the data
        minPeak = 12000;
        [t,q] = findpeaks(yData, 'NPeaks', 2, 'MinPeakProminence', minPeak);

        j = 0;
        while isempty(t)
            j = j+100;
           [t,q] = findpeaks(yData, 'NPeaks', 2, 'MinPeakProminence', minPeak-j);
        end

        % Initial guesses
        % Linewidths
        halfmax = 0.5*(max(yData)+min(yData)); % Half height
        index1 = find(yData >= halfmax, 1, 'first'); % First time above half height
        index2 = find(yData >= halfmax, 1, 'last'); % Last time above half height
        gamma = xData(index2)-xData(index1); % Approximate linewidth

        B = min(yData); % Background
        if length(q) == 1 %findpeaks() finds only a single peak
            q = [ q-1 q-1 ];
            t = [ t t ];
            x02 = xData(q(2)) + gamma/4;
            gamma1 = gamma/2;
        elseif length(q) == 2 %findpeaks() finds two peaks.
            q = [ q(1)-1 q(2)-1 ];
            x02 = xData(q(2));
            gamma1 = gamma/5;
        end

        x01 = xData(q(1));
        gamma2 = gamma1;

        % These guesses come from pluggin in the two peak points into
        % the sum of two Lorentzians, and then solving.
        A1 = (pi*gamma*(4*(x01-x02).^2+gamma.^2)*(-4*(B-t(1))*(x01-x02).^2+(t(1)-t(2))*gamma.^2))/(...
            8*(x01-x02).^2*(4*(x01-x02).^2+gamma.^2)+8*(x01-x02).^2*gamma.^2);
        A2 = -(pi*(4*(B-t(2))*(x01-x02).^2+(t(1)-t(2))*gamma.^2)*gamma*(4*(x01-x02).^2+gamma.^2))/(...
            8*(x01-x02).^2*(4*(x01-x02).^2+gamma^2)+8*(x01-x02).^2*gamma.^2);
    elseif length(guess) ~= 7
        cprintf('err', ['\nERROR: 7 parameters needed to fit, only ' num2str(length(guess)) ' supplied.\n']);
        beep; return
    else
        A1 = guess(1);
        A2 = guess(2);
        gamma1 = guess(3);
        gamma2 = guess(4);
        x01 = guess(5);
        x02 = guess(6);
        B = guess(7);
    end

    % Fit the data with a sum of two Lorentzians
    [f,gof] = fit(xData,yData,Lorentz2,'StartPoint',[A1, A2, gamma1, gamma2, x01, x02, B],...
        'Weights',max(yErr)./yErr,'Lower',[0 0 0 0 0 0 0]);
    
    % Output the guesses
    guess = [A1, A2, gamma1, gamma2, x01, x02, B];
    
end