function plots_HeNe(plotarrays)
%% plots_HeNe.m
% Make plots for scans where only the HeNe laser power was changed.

%% 
    % w0s vs. ACvalues
    figure('Name', 'w0s vs. Resonant Laser OD', 'WindowStyle', 'docked', 'numbertitle', 'off');
    
    errorbar(plotarrays.HeNeODs, plotarrays.w0s.HeNe, plotarrays.w0s.HeNeerr,...
        'r.', 'markersize', 15, 'Capsize', 0.1);
    
    xlabel('HeNe Laser OD');
    ylabel('Center Frequency (GHz)');
    
    % linewidths vs. ACvalues
    figure('Name', 'Linewidths vs. Resonant Laser OD', 'WindowStyle', 'docked', 'numbertitle', 'off');
    
    errorbar(plotarrays.HeNeODs, plotarrays.linewidths.HeNe, plotarrays.linewidths.HeNeerr,...
        'r.', 'markersize', 15,'Capsize', 0.1);
    
    xlabel('HeNe Laser OD');
    ylabel('Linewidth (GHz)');
    
    % heights-background vs. ACvalues
    figure('Name', 'Heights vs. Resonant Laser OD', 'WindowStyle', 'docked', 'numbertitle', 'off');
    
    errorbar(plotarrays.HeNeODs, plotarrays.heights.HeNe-plotarrays.B.HeNe, plotarrays.heights.HeNeerr,...
        'r.', 'markersize', 15, 'Capsize', 0.1);
    
    xlabel('HeNe Laser OD');
    ylabel('Peak Height (arb. units)');
    hold off
    
    % areas vs. ACvalues
    figure('Name', 'Areas vs. Resonant Laser OD', 'WindowStyle', 'docked', 'numbertitle', 'off');
    
    rrorbar(plotarrays.HeNeODs, plotarrays.areas.HeNe, plotarrays.areas.HeNeerr,...
        'r.', 'markersize', 15, 'Capsize', 0.1);
    
    xlabel('HeNe Laser OD');
    ylabel('Peak Area (Int*GHz)');
    
end