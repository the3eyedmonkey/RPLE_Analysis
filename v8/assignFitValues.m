function fitvalues = assignFitValues(data, fitvalues, f, numpeaks)
%% assignFitValues.m
% Output all of the fit parameters into 'fitvalues' structure
% If a single Lorentzian is used saves like [value, confint1, confint2]
% If a double Lorentzian is used saves like
% [value1, confint11, confint12, value2, confint21, confint22]
% Here
% data = data(i) from the main function
% fitvalues = fitvalues(i) from the main function
% numpeaks = numPeaks(i) from the main function

%% 
    % Record the scan name
    fitvalues.filename = data.filename;
    
    % Fit parameters and intervals
    values = coeffvalues(f);
    intervals = confint(f);
    
    if numpeaks == 1
        % Get height from fitting parameters
        [h, hint1, hint2] = heightLorentz(values, intervals);

        % Place fit parameters in the fitvalues structure
        fitvalues.heights = [ h hint1 hint2 ];
        fitvalues.areas = [ values(1) intervals(1,1) intervals(2,1) ];
        fitvalues.linewidths = [ values(2) intervals(1,2) intervals(2,2) ];
        fitvalues.w0s = [ values(3) intervals(1,3) intervals(2,3) ];
        fitvalues.B = [ values(4) intervals(1,4) intervals(2,4) ];
    elseif numpeaks == 2
        % Get height from fitting parameters
        [h1, h1int1, h1int2, h2, h2int1, h2int2] = heightsLorentz2(values, intervals);

        % Place fit parameters in the fitvalues structure
        fitvalues.heights = [ h1 h1int1 h1int2...
            h2 h2int1 h2int2 ];

        fitvalues.areas = [ values(1) intervals(1,1) intervals(2,1)...
            values(2) intervals(1,2) intervals(2,2) ];

        fitvalues.linewidths = [ values(3) intervals(1,3) intervals(2,3)...
            values(4) intervals(1,4) intervals(2,4) ];

        fitvalues.w0s = [ values(5) intervals(1,5) intervals(2,5)...
            values(6) intervals(1,6) intervals(2,6) ];

        fitvalues.B = [ values(7) intervals(1,7) intervals(2,7) ];
    end
    
end