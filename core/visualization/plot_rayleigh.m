function plot_rayleigh(U, stats, COL, config, outDir)
% Rayleigh tuning + polar plot for a single unit.

if ~isfield(stats, 'rayleigh') || isempty(stats.rayleigh)
    return;
end

RY = stats.rayleigh;
if ~isfield(RY, 'theta_deg') || ~isfield(RY, 'rate_mean')
    return;
end

theta_deg = RY.theta_deg(:);
resp      = RY.rate_mean(:);

if isempty(theta_deg) || all(~isfinite(resp))
    return;
end

% plotting flags
makePlots = true;
if isfield(config, 'plot') && isfield(config.plot, 'makePlots')
    makePlots = logical(config.plot.makePlots);
end
figVis = 'on';
if ~makePlots
    figVis = 'off';
end

pngDpi = 300;
if isfield(config, 'plot') && isfield(config.plot, 'dpi')
    pngDpi = config.plot.dpi;
end

unitDir = fullfile(outDir, sprintf('%s_%d', U.unitType, U.unitID));
if ~exist(unitDir,'dir'), mkdir(unitDir); end

% make figure with two subplots (cartesian + polar)
f = figure('Visible',figVis,'Color','w');

% --- left: tuning vs angle ---
ax1 = subplot(1,2,1, 'Parent', f);
hold(ax1,'on');

plot(ax1, theta_deg, resp, 'o-','LineWidth',1.2);
xlabel(ax1,'Hue angle (deg)');
ylabel(ax1,'Mean firing (Hz)');
grid(ax1,'on');

ttl = sprintf('Rayleigh tuning | R=%.2f, p=%.3g', RY.R, RY.p);
title(ax1, ttl);

% --- right: polar vector ---
ax2 = subplot(1,2,2, 'Parent', f);
ax2 = polaraxes('Parent', f); % put a polar axes on top of subplot 2
set(ax2, 'Position', get(ax2.Parent.Children(1),'Position')); % roughly align

hold(ax2,'on');

theta_rad = deg2rad(theta_deg);
r = resp - min(resp);
if max(r) > 0
    r = r ./ max(r);
end

polarplot(ax2, theta_rad, r, 'o-','LineWidth',1.2);

mu_rad = deg2rad(RY.mu_deg);
polarplot(ax2, [mu_rad mu_rad], [0 1], 'r-','LineWidth',2);

title(ax2, sprintf('Preferred angle = %.1f°', RY.mu_deg));

% save
fileTag = sprintf('%s_%s_%d', string(U.dateStr), U.unitType, U.unitID);
pngName = fullfile(unitDir, sprintf('Rayleigh_%s.png', fileTag));
figName = fullfile(unitDir, sprintf('Rayleigh_%s.fig', fileTag));

exportgraphics(f, pngName, 'Resolution', pngDpi);
savefig(f, figName);

if ~makePlots
    close(f);
end

end
