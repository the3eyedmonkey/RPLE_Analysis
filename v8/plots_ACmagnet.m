function plots_ACmagnet(plotarrays)
%% plots_ACmagnet.m
% Plots for scans that both the AC Stark laser and the magnet on.

%% 
    % Mean of reference scans
    refmean = mean(plotarrays.w0s.ref);
    
    % w0s
    figure('Name', 'Center Frequencies', 'WindowStyle', 'docked', 'numbertitle', 'off');
    
    % AC scans
    h1 = plot3(plotarrays.ACvalues.AC, zeros(1,length(plotarrays.ACvalues.AC)),...
        (plotarrays.w0s.AC-refmean), 'r.', 'markersize', 15);
    hold on
    % Magnet scans
    h2 = plot3(zeros(1,length(plotarrays.magvalues.mag)),plotarrays.magvalues.mag,...
        (plotarrays.w0s.mag-refmean), 'g.', 'markersize', 15);
    % ACmagnet scans
    h3 = plot3(plotarrays.ACvalues.ACmag,plotarrays.magvalues.ACmag,...
        (plotarrays.w0s.ACmag-refmean), 'm.', 'markersize', 15);
    
    legend([h1 h2 h3], 'AC', 'magnet', 'ACmagnet');
    
    xlabel('AC Stark Laser Power (mW)');
    ylabel('Magnetic Field Strength (mT)');
    zlabel('Center Frequency (GHz)');
    grid on
    hold off
    
    
    % linewidths
    figure('Name', 'Linewidths', 'WindowStyle', 'docked', 'numbertitle', 'off');
    
    % AC scans
    h1 = plot3(plotarrays.ACvalues.AC,zeros(1,length(plotarrays.ACvalues.AC)),...
        (plotarrays.linewidths.AC), 'r.', 'markersize', 15);
    hold on
    % Magnet scans
    h2 = plot3(zeros(1,length(plotarrays.magvalues.mag)), plotarrays.magvalues.mag,...
        (plotarrays.linewidths.mag), 'g.', 'markersize', 15);
    % ACmagnet scans
    h3 = plot3(plotarrays.ACvalues.ACmag, plotarrays.magvalues.ACmag,...
        (plotarrays.linewidths.ACmag), 'm.', 'markersize', 15);
    
    legend([h1 h2 h3], 'AC', 'magnet', 'ACmagnet');
    
    xlabel('AC Stark Laser Power (mW)');
    ylabel('Magnetic Field Strength (mT)');
    zlabel('Linewidth (GHz)');
    grid on
    hold off
    
    
    % heights-background
    figure('Name', 'Peak Heights', 'WindowStyle', 'docked', 'numbertitle', 'off');
    
    % AC scans
    h1 = plot3(plotarrays.ACvalues.AC, zeros(1,length(plotarrays.ACvalues.AC)),...
        (plotarrays.heights.AC-plotarrays.B.AC), 'r.', 'markersize', 15);
    hold on
    % Magnet scans
    h2 = plot3(zeros(1,length(plotarrays.magvalues.mag)), plotarrays.magvalues.mag,...
        (plotarrays.linewidths.mag-plotarrays.B.mag), 'g.', 'markersize', 15);
    % ACmagnet scans
    h3 = plot3(plotarrays.ACvalues.ACmag, plotarrays.magvalues.ACmag,...
        (plotarrays.linewidths.ACmag-plotarrays.B.ACmag), 'm.', 'markersize', 15);
    
    legend([h1 h2 h3], 'AC', 'magnet', 'ACmagnet');
    
    xlabel('AC Stark Laser Power (mW)');
    ylabel('Magnetic Field Strength (mT)');
    zlabel('Peak Height (arb. units)');
    grid on
    hold off
    
    
    % areas
    figure('Name', 'Peak Areas', 'WindowStyle', 'docked', 'numbertitle', 'off');
    
    % AC scans
    h1 = plot3(plotarrays.ACvalues.AC, zeros(1,length(plotarrays.ACvalues.AC)),...
        (plotarrays.areas.AC), 'r.', 'markersize', 15);
    hold on
    % Magnet scans
    h2 = plot3(zeros(1,length(plotarrays.magvalues.mag)), plotarrays.magvalues.mag,...
        (plotarrays.areas.mag), 'g.', 'markersize', 15);
    % ACmagnet scans
    h3 = plot3(plotarrays.ACvalues.ACmag, plotarrays.magvalues.ACmag,...
        (plotarrays.areas.ACmag), 'm.', 'markersize', 15);
    
    legend([h1 h2 h3], 'AC', 'magnet', 'ACmagnet');
    
    xlabel('AC Stark Laser Power (mW)');
    ylabel('Magnetic Field Strength (mT)');
    zlabel('Peak Area (Int*GHz)');
    grid on
    hold off
    
end