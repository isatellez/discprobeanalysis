function plot_tuning_curves(U, stats, COL, config, outDir)
% Tuning curve figure: mean rate per hue + screen colors
% loosely based on DO_TUNING_FIG block in Analyze_DiscProbePSTHs_v1figproc.m

if ~isfield(stats, 'hueMeans') || isempty(stats.hueMeans)
    return;
end

hm   = stats.hueMeans;
hues = hm.hues;
y    = hm.rate_mean;
e    = hm.rate_sem;

if isempty(hues) || all(~isfinite(y))
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

% ---------- paths: per-unit, per-session, global ----------

dateStr  = char(string(U.dateStr));
unitType = char(string(U.unitType));
unitID   = U.unitID;

% per-unit dir: units/<unit>/figures (flat)
Upaths  = get_unit_paths(config, dateStr, unitType, unitID);
unitDir = Upaths.figures;
if ~exist(unitDir,'dir')
    mkdir(unitDir);
end


% per-session dir: <date>/figs/discprobe/tuning
if nargin < 5 || isempty(outDir)
    outDir = '';
end
sessionDir = '';
if ~isempty(outDir)
    sessionDir = fullfile(outDir, 'tuning');
    if ~exist(sessionDir,'dir')
        mkdir(sessionDir);
    end
end

% global dir: output/figs/discprobe/tuning
globalDir = '';
if isfield(config,'paths') && isfield(config.paths,'globalDiscProbeFigRoot') ...
        && ~isempty(config.paths.globalDiscProbeFigRoot)
    globalDir = fullfile(config.paths.globalDiscProbeFigRoot, 'tuning');
    if ~exist(globalDir,'dir')
        mkdir(globalDir);
    end
end

% ---------- figure + plots ----------

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

    if isfield(COL, 'colsSatur') && ~isempty(COL.colsSatur)
        satIdx = size(COL.colsSatur,1); % highest saturation row
        if h >= 1 && h <= size(COL.colsSatur,2)
            rgb = squeeze(COL.colsSatur(satIdx, h, :));
            rgb = double(rgb(:)');
            if max(rgb) > 1
                rgb = rgb / 255;
            end
            rgb(~isfinite(rgb)) = 0;
            c = rgb;
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

% ---------- save in all three locations ----------

fileTag  = sprintf('%s_%s_%d', string(U.dateStr), U.unitType, U.unitID);
baseName = sprintf('TuningCurve_%s', fileTag);

% per-unit
png_unit = fullfile(unitDir, [baseName '.png']);
fig_unit = fullfile(unitDir, [baseName '.fig']);
exportgraphics(f, png_unit, 'Resolution', pngDpi);
savefig(f, fig_unit);

% per-session
if ~isempty(sessionDir)
    png_sess = fullfile(sessionDir, [baseName '.png']);
    fig_sess = fullfile(sessionDir, [baseName '.fig']);
    exportgraphics(f, png_sess, 'Resolution', pngDpi);
    savefig(f, fig_sess);
end

% global
if ~isempty(globalDir)
    png_glob = fullfile(globalDir, [baseName '.png']);
    fig_glob = fullfile(globalDir, [baseName '.fig']);
    exportgraphics(f, png_glob, 'Resolution', pngDpi);
    savefig(f, fig_glob);
end

if ~makePlots
    close(f);
end

end
