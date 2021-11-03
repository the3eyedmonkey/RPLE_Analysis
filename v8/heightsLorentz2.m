function [h1, h1int1, h1int2, h2, h2int1, h2int2] = heightsLorentz2(values, intervals)
%% heightsLorentz2.m
% Compute the hieghts from fitting parameters of two area normalized
% Lorentzians.

%% 
    %
    h1 = (2*values(1))/(pi*values(3)) + (values(2)*values(4))/(2*pi*((values(5)-values(6))^2+(0.5*values(4))^2)) + values(7);
    h1int1 = (2*intervals(1,1))/(pi*values(3)) + (values(2)*values(4))/(2*pi*((values(5)-values(6))^2+(0.5*values(4))^2)) + values(7);
    h1int2 = (2*intervals(2,1))/(pi*values(3)) + (values(2)*values(4))/(2*pi*((values(5)-values(6))^2+(0.5*values(4))^2)) + values(7);

    h2 = (values(1)*values(3))/(2*pi*((values(6)-values(5))^2+(0.5*values(3))^2)) + (2*values(2))/(pi*values(4)) + values(7);
    h2int1 = (values(1)*values(3))/(2*pi*((values(6)-values(5))^2+(0.5*values(3))^2)) + (2*intervals(1,2))/(pi*values(4)) + values(7);
    h2int2 = (values(1)*values(3))/(2*pi*((values(6)-values(5))^2+(0.5*values(3))^2)) + (2*intervals(2,2))/(pi*values(4)) + values(7);
    
end