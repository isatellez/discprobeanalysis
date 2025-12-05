function plot_fourier(U, stats, COL, config, outDir)

if ~isfield(stats, 'fourier') || isempty(stats.fourier)
    return;
end

FT = stats.fourier;
FT_over = FT.over;
FT_sat  = FT.perSat;
saturIDs = FT.saturIDs;

makePlots = true;
if isfield(config, 'plot') && isfield(config.plot, 'makePlots')
    makePlots = logical(config.plot.makePlots);
end
figVis = 'on';
if ~makePlots
    figVis = 'off';
end

unitDir = fullfile(outDir, sprintf('%s_%d', U.unitType, U.unitID));
if ~exist(unitDir,'dir'), mkdir(unitDir); end

figF = figure('Color','w','Visible',figVis,'Name','Hue Fourier');
TLf  = tiledlayout(figF,2,2,'Padding','compact','TileSpacing','compact');

if isfield(U,'exptName')
    unitLabel = sprintf('%s | Unit %d (%s) — Hue Fourier', string(U.exptName), U.unitID, U.unitType);
else
    unitLabel = sprintf('%s | Unit %d (%s) — Hue Fourier', string(U.dateStr), U.unitID, U.unitType);
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

pngDpi = 200;
if isfield(config, 'plot') && isfield(config.plot, 'dpi')
    pngDpi = config.plot.dpi;
end

fileTag = sprintf('%s_%s_%d', string(U.dateStr), U.unitType, U.unitID);
pngName = fullfile(unitDir, sprintf('FourierHue_quicklook_%s.png', fileTag));
figName = fullfile(unitDir, sprintf('FourierHue_quicklook_%s.fig', fileTag));

exportgraphics(figF, pngName, 'Resolution', pngDpi);
savefig(figF, figName);

if ~makePlots
    close(figF);
end

end
