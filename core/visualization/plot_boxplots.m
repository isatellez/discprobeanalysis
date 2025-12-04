function plot_boxplots(U, stats, COL, config, outDir)
% Box plots of early-window rates by hue × saturation
% echoes BOX PLOTS + OUTLIERS block in Analyze_DiscProbePSTHs_v1figproc.m

if U.nTrials == 0
    return;
end

win = U.winEarly;
dur = diff(win);

nTr = numel(U.spk);
rate = nan(nTr,1);
for tt = 1:nTr
    spks = U.spk{tt};
    if isempty(spks)
        rate(tt) = 0;
    else
        spks = spks(spks >= win(1) & spks < win(2));
        rate(tt) = numel(spks) ./ dur;
    end
end

hueID = U.trials.hueID;
satID = U.trials.satID;

valid = isfinite(rate) & isfinite(hueID) & isfinite(satID);
rate  = rate(valid);
hueID = hueID(valid);
satID = satID(valid);

if isempty(rate)
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

sats = unique(satID);
sats = sats(:);
nS   = numel(sats);

f = figure('Visible',figVis,'Color','w');
t = tiledlayout(f,1,nS,'TileSpacing','compact','Padding','compact');

for si = 1:nS
    s = sats(si);
    ax = nexttile(t,si);
    hold(ax,'on');

    sel = satID == s;
    if ~any(sel)
        continue;
    end

    boxplot(ax, rate(sel), hueID(sel));
    xlabel(ax,'Hue index');
    ylabel(ax,'Early rate (Hz)');
    grid(ax,'on');

    title(ax, sprintf('sat = %.2f', s));
end

fileTag = sprintf('%s_%s_%d', string(U.dateStr), U.unitType, U.unitID);
pngName = fullfile(unitDir, sprintf('Boxplots_%s.png', fileTag));
figName = fullfile(unitDir, sprintf('Boxplots_%s.fig', fileTag));

exportgraphics(f, pngName, 'Resolution', pngDpi);
savefig(f, figName);

if ~makePlots
    close(f);
end

end

function out = iff(cond, a, b)
    if cond, out = a; else, out = b; end
end
