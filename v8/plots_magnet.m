function plots_magnet(plotarrays)
%% plots_magnet.m
% Plots for scans that only have the magnet on.

%% 
    % Mean of reference scans
    refmean = mean(plotarrays.w0s.ref);
    
    % w0s vs. magvalues
    figure('Name', 'w0s(magnet only) vs. Magnetic Field Strength', 'WindowStyle', 'docked', 'numbertitle', 'off');
    hold on
    
    h1 = errorbar(plotarrays.magvalues.mag, (plotarrays.w0s.mag-refmean), plotarrays.w0s.magerr,...
        'g.', 'markersize', 15, 'Capsize', 0.1);
    h0 = hline(0, 'k:');
    
    legend([h0 h1], 'Reference', 'magnet');
    
    xlabel('Magnetic Field Strength (mT)');
    ylabel('Center Frequency (GHz)');
    hold off
    
    % linewidths vs. magvalues
    figure('Name', 'Linewidths vs. Magnetic Field Strength', 'WindowStyle', 'docked', 'numbertitle', 'off');
    hold on
    h1 = errorbar(plotarrays.magvalues.mag, plotarrays.linewidths.mag, plotarrays.linewidths.magerr,...
        'g.', 'markersize', 15, 'Capsize', 0.1);
    h0 = hline(mean(plotarrays.linewidths.ref), 'k:');
    
    legend([h0 h1], 'Reference', 'magnet');
    
    xlabel('Magnetic Field Strength (mT)');
    ylabel('Linewidth (GHz)');
    hold off
    
    % heights-background vs. magvalues
    figure('Name', 'Heights vs. Magnetic Field Strength', 'WindowStyle', 'docked', 'numbertitle', 'off');
    hold on
    h1 = errorbar(plotarrays.magvalues.mag, plotarrays.heights.mag-plotarrays.B.mag, plotarrays.heights.magerr,...
        'g.', 'markersize', 15, 'Capsize', 0.1);
    h0 = hline(mean(plotarrays.heights.ref-plotarrays.B.ref), 'k:');
    
    legend([h0 h1], 'Reference', 'magnet');
    
    xlabel('Magnetic Field Strength (mT)');
    ylabel('Peak Height (GHz)');
    hold off
    
    % areas vs. magvalues
    figure('Name', 'Areas vs. Magnetic Field Strength', 'WindowStyle', 'docked', 'numbertitle', 'off');
    hold on
    h1 = errorbar(plotarrays.magvalues.mag, plotarrays.areas.mag, plotarrays.areas.magerr,...
        'g.', 'markersize', 15, 'Capsize', 0.1);
    h0 = hline(mean(plotarrays.areas.ref),'k:');
    
    legend([h0 h1], 'Reference', 'magnet');
    
    xlabel('Magnetic Field Strength (mT)');
    ylabel('Peak Area (Int*GHz)');
    hold off
    
end