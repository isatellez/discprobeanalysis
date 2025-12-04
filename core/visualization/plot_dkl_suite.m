function plot_dkl_suite(U, stats, COL, config, outDir)
% DKL visualization: L–M scatter colored by rate + polar direction map
% loosely mirrors the DKL suite in Analyze_DiscProbePSTHs_v1figproc.m

if ~all(ismember({'dklL','dklM'}, U.trials.Properties.VariableNames))
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

dL = U.trials.dklL;
dM = U.trials.dklM;

valid = isfinite(rate) & isfinite(dL) & isfinite(dM);
if ~any(valid)
    return;
end

dL = dL(valid);
dM = dM(valid);
rate = rate(valid);

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
t = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');

% L–M plane scatter
ax1 = nexttile(t,1);
scatter(ax1, dL, dM, 25, rate, 'filled');
axis(ax1,'equal');
grid(ax1,'on');
xlabel(ax1,'DKL L');
ylabel(ax1,'DKL M');
cb = colorbar(ax1);
cb.Label.String = 'Early rate (Hz)';
title(ax1, 'DKL L–M scatter');

% polar version of same info
ax2 = polaraxes(f);
ax2.Layout.Tile = 2;
hold(ax2,'on');

theta = atan2(dM, dL);
r = sqrt(dL.^2 + dM.^2);

if max(r) > 0
    r_norm = r ./ max(r);
else
    r_norm = r;
end

polarscatter(ax2, theta, r_norm, 25, rate, 'filled');
title(ax2, 'DKL direction (normalized radius)');

fileTag = sprintf('%s_%s_%d', string(U.dateStr), U.unitType, U.unitID);
pngName = fullfile(unitDir, sprintf('DKL_suite_%s.png', fileTag));
figName = fullfile(unitDir, sprintf('DKL_suite_%s.fig', fileTag));

exportgraphics(f, pngName, 'Resolution', pngDpi);
savefig(f, figName);

if ~makePlots
    close(f);
end

end

function out = iff(cond, a, b)
    if cond, out = a; else, out = b; end
end
