function plot_fourier(U, stats, COL, config, outDir)

if ~isfield(stats, 'fourier') || isempty(stats.fourier)
    return;
end

FT       = stats.fourier;
FT_over  = FT.over;
FT_sat   = FT.perSat;
saturIDs = FT.saturIDs;

% figure visibility
makePlots = true;
if isfield(config, 'plot') && isfield(config.plot, 'makePlots')
    makePlots = logical(config.plot.makePlots);
end
figVis = 'on';
if ~makePlots
    figVis = 'off';
end

% ---------- paths: per-unit, per-session, global ----------

dateStr  = char(string(U.dateStr));
unitType = char(string(U.unitType));
unitID   = U.unitID;


% per-unit dir: units/<unit>/figures/discprobe/fourier
Upaths  = get_unit_paths(config, dateStr, unitType, unitID);
unitDir = Upaths.figures;
if ~exist(unitDir,'dir')
    mkdir(unitDir);
end


% per-session dir: <date>/figs/discprobe/fourier
if nargin < 5 || isempty(outDir)
    outDir = '';
end
sessionDir = '';
if ~isempty(outDir)
    sessionDir = fullfile(outDir, 'fourier');
    if ~exist(sessionDir,'dir')
        mkdir(sessionDir);
    end
end

% global dir: output/figs/discprobe/fourier
globalDir = '';
if isfield(config,'paths') && isfield(config.paths,'globalDiscProbeFigRoot') ...
        && ~isempty(config.paths.globalDiscProbeFigRoot)
    globalDir = fullfile(config.paths.globalDiscProbeFigRoot, 'fourier');
    if ~exist(globalDir,'dir')
        mkdir(globalDir);
    end
end

% ---------- plotting ----------

figF = figure('Color','w','Visible',figVis,'Name','Hue Fourier');
TLf  = tiledlayout(figF,2,2,'Padding','compact','TileSpacing','compact');

if isfield(U,'exptName')
    unitLabel = sprintf('%s | Unit %d (%s) — Hue Fourier', ...
        string(U.exptName), U.unitID, U.unitType);
else
    unitLabel = sprintf('%s | Unit %d (%s) — Hue Fourier', ...
        string(U.dateStr), U.unitID, U.unitType);
end
sgtitle(figF, unitLabel, 'Interpreter','none');

% amplitude spectrum (overall)
ax = nexttile(TLf);
bar(ax, FT_over.harmonics, FT_over.amp);
grid(ax,'on');
xlabel(ax,'Harmonic k');
ylabel(ax,'Amplitude');
title(ax,'Amplitude spectrum');

% variance explained (overall, no-DC fraction)
ax = nexttile(TLf);
bar(ax, FT_over.harmonics, FT_over.varExp);
grid(ax,'on');
xlabel(ax,'Harmonic k');
ylabel(ax,'Variance explained (noDC frac)');
title(ax,'Variance explained');

% first harmonic amplitude per sat
nS = numel(FT_sat);
ax = nexttile(TLf);
a1 = NaN(1, nS);
for si = 1:nS
    if isfield(FT_sat(si),'amp') && numel(FT_sat(si).amp) >= 1
        a1(si) = FT_sat(si).amp(1);
    end
end
bar(ax, 1:nS, a1);
set(ax,'XTick',1:nS,'XTickLabel',compose('sat=%.2f',saturIDs));
ylabel(ax,'Amp k=1');
title(ax,'First-harmonic amplitude (per sat)');
grid(ax,'on');

% first harmonic phase per sat
ax = nexttile(TLf);
ph1 = NaN(1, nS);
for si = 1:nS
    if isfield(FT_sat(si),'phase_deg') && numel(FT_sat(si).phase_deg) >= 1
        ph1(si) = FT_sat(si).phase_deg(1);
    end
end
bar(ax, 1:nS, ph1);
set(ax,'XTick',1:nS,'XTickLabel',compose('sat=%.2f',saturIDs));
ylabel(ax,'Phase k=1 (deg)');
title(ax,'First-harmonic phase (per sat)');
grid(ax,'on');

% ---------- save to all three locations ----------

pngDpi = 200;
if isfield(config, 'plot') && isfield(config.plot, 'dpi')
    pngDpi = config.plot.dpi;
end

fileTag = sprintf('%s_%s_%d', string(U.dateStr), U.unitType, U.unitID);

% per-unit
png_unit = fullfile(unitDir, sprintf('FourierHue_quicklook_%s.png', fileTag));
fig_unit = fullfile(unitDir, sprintf('FourierHue_quicklook_%s.fig', fileTag));
exportgraphics(figF, png_unit, 'Resolution', pngDpi);
savefig(figF, fig_unit);

% per-session
if ~isempty(sessionDir)
    png_sess = fullfile(sessionDir, sprintf('FourierHue_quicklook_%s.png', fileTag));
    fig_sess = fullfile(sessionDir, sprintf('FourierHue_quicklook_%s.fig', fileTag));
    exportgraphics(figF, png_sess, 'Resolution', pngDpi);
    savefig(figF, fig_sess);
end

% global
if ~isempty(globalDir)
    png_glob = fullfile(globalDir, sprintf('FourierHue_quicklook_%s.png', fileTag));
    fig_glob = fullfile(globalDir, sprintf('FourierHue_quicklook_%s.fig', fileTag));
    exportgraphics(figF, png_glob, 'Resolution', pngDpi);
    savefig(figF, fig_glob);
end

if ~makePlots
    close(figF);
end

end
