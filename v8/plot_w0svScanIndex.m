function plot_w0svScanIndex(plotarrays)
%% plot_w0svScanIndex.m
% Name says it all.

%% 
    % Mean of reference scans
    refmean = mean(plotarrays.w0s.ref);
    
    % Make a plot of w0s, with error bars for reference
    figure('Name', 'w0s vs. Scan Index', 'WindowStyle', 'docked', 'numbertitle', 'off');
    hold on
    
    errorbar(plotarrays.indices.ref, (plotarrays.w0s.ref-refmean), plotarrays.w0s.referr,...
        'b.', 'markersize', 15, 'Capsize', 0.1)
    errorbar(plotarrays.indices.AC, (plotarrays.w0s.AC-refmean), plotarrays.w0s.ACerr,...
        'r.', 'markersize', 15, 'Capsize', 0.1)
    errorbar(plotarrays.indices.mag, (plotarrays.w0s.mag-refmean), plotarrays.w0s.magerr,...
        'g.', 'markersize', 15, 'Capsize', 0.1)
    errorbar(plotarrays.indices.ACmag, (plotarrays.w0s.ACmag-refmean), plotarrays.w0s.ACmagerr,...
        'm.', 'markersize', 15, 'Capsize', 0.1)
    hline(0, 'k:')
    
    legend('Reference', 'AC', 'magnet', 'ACmagnet');
    
    xlabel('Scan Index');
    ylabel('w0s (GHz)');
    hold off
    
end