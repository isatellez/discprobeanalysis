function H = plot_tuning_population_histograms(pop, config)
% Histograms of preferred direction and tuning width


set(0,'DefaultFigureVisible','off');

    phi  = pop.phi0_deg(pop.good);
    fwhm = pop.fwhm_deg(pop.good);

    H.fig = figure('Color','w');

    subplot(1,2,1);
    edges = 0:30:360;
    histogram(phi, edges);
    xlabel('Preferred direction (deg)');
    ylabel('Count');
    title('Tuning peak directions');

    subplot(1,2,2);
    edgesW = 0:15:180;
    histogram(fwhm, edgesW);
    xlabel('FWHM (deg)');
    ylabel('Count');
    title('Tuning widths');

    % optional: draw vertical line for cosine width reference if you want
end
