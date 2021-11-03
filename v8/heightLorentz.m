function [h, hint1, hint2] = heightLorentz(values, intervals)
%% heightLorentz.m
% Compute the hieght from fitting parameters of area normalized
% Lorentzian.

%% 
    %
    h = (2*values(1))/(pi*values(2)) + values(4);
    hint1 = (2*intervals(1,1))/(pi*values(2)) + values(4);
    hint2 = (2*intervals(2,1))/(pi*values(2)) + values(4);
    
end