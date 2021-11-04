function plots_res(plotarrays)
%% plots_res.m
% Make plots for scans where only the resonant laser power was changed.

%% 
    % w0s vs. ACvalues
    figure('Name', 'w0s vs. Resonant Laser OD', 'WindowStyle', 'docked', 'numbertitle', 'off');
    
    errorbar(plotarrays.resODs, plotarrays.w0s.res, plotarrays.w0s.reserr,...
        'r.', 'markersize', 15, 'Capsize', 0.1);
    
    xlabel('Resonant Laser OD');
    ylabel('Center Frequency (GHz)');
    
    % linewidths vs. ACvalues
    figure('Name', 'Linewidths vs. Resonant Laser OD', 'WindowStyle', 'docked', 'numbertitle', 'off');
    
    errorbar(plotarrays.resODs, plotarrays.linewidths.res, plotarrays.linewidths.reserr,...
        'r.', 'markersize', 15,'Capsize', 0.1);
    
    xlabel('Resonant Laser OD');
    ylabel('Linewidth (GHz)');
    
    % heights-background vs. ACvalues
    figure('Name', 'Heights vs. Resonant Laser OD', 'WindowStyle', 'docked', 'numbertitle', 'off');
    
    errorbar(plotarrays.resODs, plotarrays.heights.res-plotarrays.B.res, plotarrays.heights.reserr,...
        'r.', 'markersize', 15, 'Capsize', 0.1);
    
    xlabel('Resonant Laser OD');
    ylabel('Peak Height (arb. units)');
    hold off
    
    % areas vs. ACvalues
    figure('Name', 'Areas vs. Resonant Laser OD', 'WindowStyle', 'docked', 'numbertitle', 'off');
    
    errorbar(plotarrays.resODs, plotarrays.areas.res, plotarrays.areas.reserr,...
        'r.', 'markersize', 15, 'Capsize', 0.1);
    
    xlabel('Resonant Laser OD');
    ylabel('Peak Area (Int*GHz)');
    
end