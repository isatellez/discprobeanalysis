function plot_spike_accounting(U, stats, COL, config, outDir)
% Spike accounting figure: spikes/trials per hue, grouped by sat, cosine fit
% based on SPIKE ACCOUNTING section in Analyze_DiscProbePSTHs_v1figproc.m

if U.nTrials == 0
    return;
end

win = U.winEarly;
dur = diff(win);

nTr = numel(U.spk);
spkCount = nan(nTr,1);
rate     = nan(nTr,1);

for tt = 1:nTr
    spks = U.spk{tt};
    if isempty(spks)
        spkCount(tt) = 0;
        rate(tt)     = 0;
    else
        spks = spks(spks >= win(1) & spks < win(2));
        spkCount(tt) = numel(spks);
        rate(tt)     = spkCount(tt) ./ dur;
    end
end

hueID = U.trials.hueID;
satID = U.trials.satID;

valid = isfinite(hueID) & isfinite(satID) & isfinite(rate);
hueID = hueID(valid);
satID = satID(valid);
spkCount = spkCount(valid);
rate     = rate(valid);

if isempty(hueID)
    return;
end

hues = unique(hueID);
hues = hues(:);
nHue = max(hues);

% spikes per hue
spikesPerHue = accumarray(hueID, spkCount, [nHue 1], @sum, 0);

% trials per hue and sat
sats = unique(satID);
sats = sats(:);
nS   = numel(sats);

trialCounts = accumarray([hueID, satID_index(satID, sats)], 1, [nHue nS], @sum, 0);

% mean rate per hue from stats
if isfield(stats, 'hueMeans') && isfield(stats.hueMeans, 'rate_mean')
    R_h = nan(nHue,1);
    hm = stats.hueMeans;
    R_h(hm.hues) = hm.rate_mean;
else
    R_h = accumarray(hueID, rate, [nHue 1], @mean, NaN);
end

R_norm = R_h;
R_norm = R_norm - min(R_norm,[],'omitnan');
den = max(R_norm,[],'omitnan');
if den > 0
    R_norm = R_norm ./ den;
end

% cosine fit if available
hasCos = isfield(stats, 'cosine') && isfield(stats.cosine, 'y_hat');
if hasCos
    y_hat = stats.cosine.y_hat;
    if numel(y_hat) ~= nHue
        hasCos = false;
    end
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
t = tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');

% spikes per hue
ax1 = nexttile(t,1);
bar(ax1, 1:nHue, spikesPerHue);
xlabel(ax1,'Hue index');
ylabel(ax1,'Spikes (early window)');
title(ax1,'Total spikes per hue');
grid(ax1,'on');

% trials per hue × sat
ax2 = nexttile(t,2);
hb = bar(ax2, trialCounts, 'grouped');
xlabel(ax2,'Hue index');
ylabel(ax2,'Trials');
title(ax2,'Trials per hue (by sat)');
grid(ax2,'on');

legend(ax2, compose('sat=%.2f', sats), 'Location','northeastoutside');

% normalized mean + cosine fit
ax3 = nexttile(t,3);
b3 = bar(ax3, 1:nHue, R_norm);
b3.FaceAlpha = 0.7;
hold(ax3,'on');

if hasCos
    hx = 1:nHue;
    yhat_norm = y_hat;
    yhat_norm = yhat_norm - min(yhat_norm,[],'omitnan');
    d2 = max(yhat_norm,[],'omitnan');
    if d2 > 0
        yhat_norm = yhat_norm ./ d2;
    end
    plot(ax3, hx, yhat_norm, 'k--','LineWidth',2);
end

xlabel(ax3,'Hue index');
ylabel(ax3,'Normalized mean rate');
title(ax3,'Normalized tuning ± cosine fit');
grid(ax3,'on');

% bottom row: throw summary text
ax4 = nexttile(t, [1 3]);
axis(ax4,'off');

txt = {
    sprintf('Date: %s', string(U.dateStr))
    sprintf('Unit %d (%s)', U.unitID, U.unitType)
    sprintf('nTrials = %d', numel(rate))
    sprintf('Total spikes (early window) = %d', sum(spkCount))
    };

if hasCos && isfield(stats.cosine, 'amp') && isfield(stats.cosine, 'pref_deg')
    txt{end+1} = sprintf('Cosine amp = %.3f', stats.cosine.amp);
    txt{end+1} = sprintf('Pref angle = %.1f°', stats.cosine.pref_deg);
end

if isfield(stats, 'permutation') && isfield(stats.permutation,'p')
    txt{end+1} = sprintf('Permutation p (tuning structure) = %.3g', stats.permutation.p);
end

text(ax4, 0.01, 0.9, txt, 'Units','normalized','VerticalAlignment','top');

title(ax4,'Spike accounting summary','FontWeight','bold');

fileTag = sprintf('%s_%s_%d', string(U.dateStr), U.unitType, U.unitID);
pngName = fullfile(unitDir, sprintf('SpikeAccounting_%s.png', fileTag));
figName = fullfile(unitDir, sprintf('SpikeAccounting_%s.fig', fileTag));

exportgraphics(f, pngName, 'Resolution', pngDpi);
savefig(f, figName);

if ~makePlots
    close(f);
end

end

function idx = satID_index(satID, sats)
    [~, idx] = ismember(satID, sats);
    idx(idx == 0) = 1;  % should never really happen if sats built from satID
end

function out = iff(cond, a, b)
    if cond, out = a; else, out = b; end
end
