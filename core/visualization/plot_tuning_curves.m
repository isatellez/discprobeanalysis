function plot_tuning_curves(U, stats, COL, config, outDir)
% Tuning curve figure: mean rate per hue + screen colors
% loosely based on DO_TUNING_FIG block in Analyze_DiscProbePSTHs_v1figproc.m

if ~isfield(stats, 'hueMeans') || isempty(stats.hueMeans)
    return;
end

hm = stats.hueMeans;
hues = hm.hues;
y    = hm.rate_mean;
e    = hm.rate_sem;

if isempty(hues) || all(~isfinite(y))
    return;
end

makePlots = true;
if isfield(config, 'plot') && isfield(config.plot, 'makePlots')
    makePlots = logical(config.plot.makePlots);
end
figVis = iff(makePlots,'on','off');

pngDpi = 300;
if isfield(config, 'plot') && isfield(config.plot, 'dpi')
    pngDpi = config.plot.dpi;
end

unitDir = fullfile(outDir, sprintf('%s_%d', U.unitType, U.unitID));
if ~exist(unitDir,'dir'), mkdir(unitDir); end

f = figure('Visible',figVis,'Color','w');
t = tiledlayout(f, 2, 1, 'TileSpacing','compact', 'Padding','compact');

% tuning panel
ax1 = nexttile(t,1);
hold(ax1,'on');

if all(isfinite(e))
    errorbar(ax1, hues, y, e, 'o-', 'LineWidth',1.2);
else
    plot(ax1, hues, y, 'o-', 'LineWidth',1.2);
end

xlabel(ax1,'Hue index');
ylabel(ax1,'Mean firing (Hz)');
set(ax1,'XTick',hues);
grid(ax1,'on');

title(ax1, sprintf('Color tuning | %s unit %d (%s)', ...
    string(U.dateStr), U.unitID, U.unitType), ...
    'Interpreter','none');

% color swatches
ax2 = nexttile(t,2);
hold(ax2,'on');

nH = numel(hues);
for i = 1:nH
    h = hues(i);
    c = [0.5 0.5 0.5];
    if isfield(COL, 'colsSatur')
        satIdx = size(COL.colsSatur,1); % highest saturation row
        if h >= 1 && h <= size(COL.colsSatur,2)
            rgb = squeeze(COL.colsSatur(satIdx, h, :));
            c = double(rgb(:)') ./ 255;
        end
    end
    patch(ax2, [i-0.4 i+0.4 i+0.4 i-0.4], [0 0 1 1], c, ...
        'EdgeColor','none');
end

xlim(ax2,[min(hues)-0.5 max(hues)+0.5]);
ylim(ax2,[0 1]);
set(ax2,'XTick',hues,'YTick',[]);
xlabel(ax2,'Hue index');
title(ax2, 'Screen colors (max sat)');

fileTag = sprintf('%s_%s_%d', string(U.dateStr), U.unitType, U.unitID);
pngName = fullfile(unitDir, sprintf('TuningCurve_%s.png', fileTag));
figName = fullfile(unitDir, sprintf('TuningCurve_%s.fig', fileTag));

exportgraphics(f, pngName, 'Resolution', pngDpi);
savefig(f, figName);

if ~makePlots
    close(f);
end

end

function out = iff(cond, a, b)
    if cond, out = a; else, out = b; end
end
