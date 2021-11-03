function plots_AC(plotarrays)
%% plots_AC.m
% Plots for scans that only have the AC Stark laser on.

%% 
    % Mean of reference scans
    refmean = mean(plotarrays.w0s.ref);
    
    % w0s vs. ACvalues
    figure('Name', 'w0s(AC only) vs. AC Stark Laser Power', 'WindowStyle', 'docked', 'numbertitle', 'off');
    hold on
    
    h1 = errorbar(plotarrays.ACvalues.AC, (plotarrays.w0s.AC-refmean), plotarrays.w0s.ACerr,...
        'r.', 'markersize', 15, 'Capsize', 0.1);
    h0 = hline(0, 'k:');
    
    legend([h0 h1], 'Reference', 'AC');
    
    xlabel('AC Stark Laser Power (mW)');
    ylabel('Center Frequency (GHz)');
    hold off
    
    % linewidths vs. ACvalues
    figure('Name', 'Linewidths vs. AC Stark Laser Power', 'WindowStyle', 'docked', 'numbertitle', 'off');
    hold on
    h1 = errorbar(plotarrays.ACvalues.AC, plotarrays.linewidths.AC, plotarrays.linewidths.ACerr,...
        'r.', 'markersize', 15,'Capsize', 0.1);
    h0 = hline(mean(plotarrays.linewidths.ref), 'k:');
    
    legend([h0 h1], 'Reference', 'AC');
    
    xlabel('AC Stark Laser Power (mW)');
    ylabel('Linewidth (GHz)');
    hold off
    
    % heights-background vs. ACvalues
    figure('Name', 'Heights vs. AC Stark Laser Power', 'WindowStyle', 'docked', 'numbertitle', 'off');
    hold on
    h1 = errorbar(plotarrays.ACvalues.AC, plotarrays.heights.AC-plotarrays.B.AC, plotarrays.heights.ACerr,...
        'r.', 'markersize', 15, 'Capsize', 0.1);
    h0 = hline(mean(plotarrays.heights.ref-plotarrays.B.ref), 'k:');
    
    legend([h0 h1], 'Reference', 'AC');
    
    xlabel('AC Stark Laser Power (mW)');
    ylabel('Peak Height (arb. units)');
    hold off
    
    % areas vs. ACvalues
    figure('Name', 'Areas vs. AC Stark Laser Power', 'WindowStyle', 'docked', 'numbertitle', 'off');
    hold on
    h1 = errorbar(plotarrays.ACvalues.AC, plotarrays.areas.AC, plotarrays.areas.ACerr,...
        'r.', 'markersize', 15, 'Capsize', 0.1);
    h0 = hline(mean(plotarrays.areas.ref), 'k:');
    
    legend([h0 h1], 'Reference', 'AC');
    
    xlabel('AC Stark Laser Power (mW)');
    ylabel('Peak Area (Int*GHz)');
    hold off
    
end