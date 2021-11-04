function plots_resHeNe(plotarrays)
%% plots_resHeNe.m
% Make plots for when both the resonant and HeNe laser ODs are being
% changed. Right now this is just a 3D plot for each of the fit parameters.

%% 
    % w0s
    figure('Name', 'Center Frequencies', 'WindowStyle', 'docked', 'numbertitle', 'off');
    
    plot3(plotarrays.resODs,plotarrays.HeNeODs,...
        plotarrays.w0s.resHeNe, 'r.', 'markersize', 15);
    
    xlabel('Resonant Laser OD');
    ylabel('HeNe Laser OD');
    zlabel('Center Frequency (GHz)');
    grid on
    
    % linewidths
    figure('Name', 'Linewidths', 'WindowStyle', 'docked', 'numbertitle', 'off');
    
    plot3(plotarrays.resODs,plotarrays.HeNeODs,...
        plotarrays.linewidths.resHeNe, 'r.', 'markersize', 15);
    
    xlabel('Resonant Laser OD');
    ylabel('HeNe Laser OD');
    zlabel('Linewidth (GHz)');
    grid on
    
    % heights-background
    figure('Name', 'Peak Heights', 'WindowStyle', 'docked', 'numbertitle', 'off');
    
    plot3(plotarrays.resODs,plotarrays.HeNeODs,...
        plotarrays.heights.resHeNe-plotarrays.B.resHeNe, 'r.', 'markersize', 15);
    
    xlabel('Resonant Laser OD');
    ylabel('HeNe Laser OD');
    zlabel('Peak Height (arb. units)');
    grid on
    
    % areas
    figure('Name', 'Peak Areas', 'WindowStyle', 'docked', 'numbertitle', 'off');
    
    plot3(plotarrays.resODs,plotarrays.HeNeODs,...
        plotarrays.areas.resHeNe, 'r.', 'markersize', 15);
    
    xlabel('Resonant Laser OD');
    ylabel('HeNe Laser OD');
    zlabel('Peak Area (Int*GHz)');
    grid on
    
end