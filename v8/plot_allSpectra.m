function plot_allSpectra(data)
%% plot_allSpectra.m
% Name says it all.

%% 
    % Total number of scans
    N = length(data);
    
    figure('Name', 'All Scans', 'WindowStyle', 'docked', 'numbertitle', 'off');
    hold on
    
    cmap = jet(N);
    
    % Plot raw data
    for i = 1:N
        plot(data(i).x,data(i).y,'Color',cmap(i,:))
    end
    
    % Label spectra with filenames
    legend(data(1:end).filename)
    
    xlabel('Frequency (GHz)');
    ylabel('Intensity (arb. units)');
    hold off