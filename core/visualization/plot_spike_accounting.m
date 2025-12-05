function plot_spike_accounting(U, stats, COL, config, outDir)

if U.nTrials == 0
    return;
end

% plotting flags
makePlots = true;
if isfield(config,'plot') && isfield(config.plot,'makePlots')
    makePlots = logical(config.plot.makePlots);
end
figVis = 'on';
if ~makePlots
    figVis = 'off';
end

pngDpi = 300;
if isfield(config,'plot') && isfield(config.plot,'dpi')
    pngDpi = config.plot.dpi;
end

unitDir = fullfile(outDir, sprintf('%s_%d', U.unitType, U.unitID));
if ~exist(unitDir,'dir'), mkdir(unitDir); end

% -------------------------------------------------------------------------
% trial-level early window firing (same as boxplots / outliers)
% -------------------------------------------------------------------------
win = U.winEarly;
dur = diff(win);

nTr = numel(U.spk);
rateEarly = nan(nTr,1);
for tt = 1:nTr
    spks = U.spk{tt};
    if isempty(spks)
        rateEarly(tt) = 0;
    else
        spks = spks(spks >= win(1) & spks < win(2));
        rateEarly(tt) = numel(spks) / dur;
    end
end

% normalize 0–1 like spk_early_norm in the legacy script
r = rateEarly;
r = r - min(r);
denom = max(r);
if denom <= 0
    spk_early_norm = zeros(size(r));
else
    spk_early_norm = r ./ denom;
end

hueID = U.trials.hueID;
satID = U.trials.satID;

% basic hue / sat info
if isfield(COL,'nHue')
    nHue = COL.nHue;
else
    nHue = max(hueID(~isnan(hueID)));
end

ZERO_HUE = 1;
if isfield(config,'space') && isfield(config.space,'zeroHue')
    ZERO_HUE = config.space.zeroHue;
end

saturIDs = unique(satID(~isnan(satID)));
saturIDs = saturIDs(:).';
nS       = numel(saturIDs);

[~, satIdx_u] = ismember(satID, saturIDs);
hueIdx_u      = hueID;

% -------------------------------------------------------------------------
% mean normalized rate per hue and per hue×sat (R_h, R_hs) + N_hs
% -------------------------------------------------------------------------
R_h  = nan(nHue,1);
N_h  = zeros(nHue,1);

for h = 1:nHue
    sel = (hueIdx_u == h);
    vals = spk_early_norm(sel);
    R_h(h) = mean(vals, 'omitnan');
    N_h(h) = sum(sel);
end

R_hs = nan(nHue, nS);
N_hs = zeros(nHue, nS);

for si = 1:nS
    sVal = saturIDs(si);
    selS = (satIdx_u == si);
    for h = 1:nHue
        sel = selS & (hueIdx_u == h);
        vals = spk_early_norm(sel);
        R_hs(h,si) = mean(vals, 'omitnan');
        N_hs(h,si) = sum(sel);
    end
end

% clean NaNs for plotting / cosine fits (same as legacy)
R_h(~isfinite(R_h)) = 0;
for si = 1:size(R_hs,2)
    col = R_hs(:,si);
    col(~isfinite(col)) = 0;
    R_hs(:,si) = col;
end

% color table per sat×hue (cols_satur2)
cols_satur2 = [];
if isfield(COL,'colsSatur') && ~isempty(COL.colsSatur)
    sz = size(COL.colsSatur);                  % nSat x nHue x 3
    cols_satur2 = reshape(double(COL.colsSatur), sz(1), sz(2), 3);
end

% -------------------------------------------------------------------------
% figure layout copied from DO_SPIKE_ACCOUNTING
% -------------------------------------------------------------------------
rootVis0 = get(groot,'DefaultFigureVisible');
set(groot,'DefaultFigureVisible','off');
cleanupFigVis = onCleanup(@() set(groot,'DefaultFigureVisible', rootVis0));

figAcc = figure('Name','Spike accounting','Color','w','Visible',figVis);
figAcc.Units    = 'normalized';
figAcc.Position = [0.1 0.1 0.8 0.8];

LM = 0.08; RM = 0.035; GAPX = 0.065;
Wax = (1 - LM - RM - 2*GAPX) / 3;
Xs  = [LM, LM+Wax+GAPX, LM+2*(Wax+GAPX)];
footerH      = 0.045;
yFooter      = 0.018;
labH         = 0.030;
gapAftFooter = 0.018;
gapRows      = 0.085;
axH          = 0.295;

yBot = yFooter + footerH + gapAftFooter + labH;
yTop = yBot + axH + gapRows + labH;

ax1 = axes('Parent',figAcc,'Units','normalized','Position',[Xs(1) yTop Wax axH]);
ax2 = axes('Parent',figAcc,'Units','normalized','Position',[Xs(2) yTop Wax axH]);
ax3 = axes('Parent',figAcc,'Units','normalized','Position',[Xs(3) yTop Wax axH]);
ax4 = axes('Parent',figAcc,'Units','normalized','Position',[Xs(1) yBot Wax axH]);
ax5 = axes('Parent',figAcc,'Units','normalized','Position',[Xs(2) yBot Wax axH]);
ax6 = axes('Parent',figAcc,'Units','normalized','Position',[Xs(3) yBot Wax axH]);

axSat = [ax3, ax4, ax5, ax6];   % overall + first 3 sats, like legacy

% -------------------------------------------------------------------------
% top-left: grouped mean normalized rate by hue×sat
% -------------------------------------------------------------------------
axes(ax1); cla(ax1);
bh = bar(ax1, R_hs, 'grouped'); hold(ax1,'on');
if ~isempty(cols_satur2)
    nSatCols = size(cols_satur2,1);
    for s = 1:min(nS, numel(bh))
        satRow = min(s, nSatCols);
        c = squeeze(cols_satur2(satRow,1:nHue,:));
        c = double(c);
        if max(c(:)) > 1, c = c/255; end
        bh(s).FaceColor = 'flat';
        bh(s).CData     = c;
    end
end
xlabel(ax1,'Hue');
ylabel(ax1,'Normalized mean spike rate (0–1)');
lg1 = legend(ax1, compose('sat=%.2f',saturIDs), ...
    'Location','northeastoutside','AutoUpdate','off');
if ~makePlots
    set(lg1,'Visible','off');
end
title(ax1,'Normalized mean spike rate per hue (grouped)');
grid(ax1,'on');
xlim(ax1,[0.5, nHue+0.5]);

% -------------------------------------------------------------------------
% top-mid: grouped trial counts per hue×sat
% -------------------------------------------------------------------------
axes(ax2); cla(ax2);
bt = bar(ax2, N_hs, 'grouped'); hold(ax2,'on');
if ~isempty(cols_satur2)
    nSatCols = size(cols_satur2,1);
    for s = 1:min(nS, numel(bt))
        satRow = min(s, nSatCols);
        c = squeeze(cols_satur2(satRow,1:nHue,:));
        c = double(c);
        if max(c(:)) > 1, c = c/255; end
        bt(s).FaceColor = 'flat';
        bt(s).CData     = c;
    end
end
xlabel(ax2,'Hue');
ylabel(ax2,'Trials');
lg2 = legend(ax2, compose('sat=%.2f',saturIDs), ...
    'Location','northeastoutside','AutoUpdate','off');
if ~makePlots
    set(lg2,'Visible','off');
end
title(ax2,'Trials per hue (grouped)');
grid(ax2,'on');
xlim(ax2,[0.5, nHue+0.5]);

% -------------------------------------------------------------------------
% cosine fit helper (same model as legacy fit_sine_hues)
% -------------------------------------------------------------------------
[aT, ampT, phT, RsqT] = fit_sine_hues_local(R_h);

stepDeg  = 360 / nHue;
shiftDeg = stepDeg * (ZERO_HUE - 1);
peakT    = mod(-phT - shiftDeg, 360);

% -------------------------------------------------------------------------
% top-right: overall normalized mean + cosine fit
% -------------------------------------------------------------------------
axes(ax3); cla(ax3);
b3 = bar(ax3, 1:nHue, R_h); hold(ax3,'on');
b3.FaceColor = 'flat';
if ~isempty(cols_satur2)
    satRow = size(cols_satur2,1);   % highest sat row, like legacy
    c = squeeze(cols_satur2(satRow,1:nHue,:));
    c = double(c);
    if max(c(:)) > 1, c = c/255; end
    b3.CData = c;
end
xlabel(ax3,'Hue');
ylabel(ax3,'Normalized mean spike rate (0–1)');
title(ax3,'Normalized mean spike rate (overall)');
grid(ax3,'on');
xlim(ax3,[0.5, nHue+0.5]);

hx = linspace(1, nHue, 300);
plot(ax3, hx, aT + ampT*cos(2*pi*(hx-1)/nHue + deg2rad(phT)), ...
    'k--', 'LineWidth', 2);

txtPeak = sprintf('peak ≈ %.0f°', peakT);
yl = ylim(ax3);
text(ax3, 0.02, 0.95, txtPeak, ...
    'Units','normalized','HorizontalAlignment','left', ...
    'VerticalAlignment','top');

% -------------------------------------------------------------------------
% bottom row: per-sat mean + cosine fit (first 3 sats)
% -------------------------------------------------------------------------
nPlotSat = min(nS, 3);
for si = 1:nPlotSat
    ax  = axSat(si+1);
    ysi = R_hs(:,si);

    axes(ax); cla(ax);
    b = bar(ax, 1:nHue, ysi); hold(ax,'on');

    b.FaceColor = 'flat';
    if ~isempty(cols_satur2)
        nSatCols = size(cols_satur2,1);
        satRow   = min(si, nSatCols);
        c = squeeze(cols_satur2(satRow,1:nHue,:));
        c = double(c);
        if max(c(:)) > 1, c = c/255; end
        b.CData = c;
    end

    xlabel(ax,'Hue');
    ylabel(ax,'Normalized mean spike rate (0–1)');
    title(ax, sprintf('sat=%.2f', saturIDs(si)));
    grid(ax,'on');
    xlim(ax,[0.5, nHue+0.5]);

    [ai,ampi,phi,Rsqi] = fit_sine_hues_local(ysi);

    hx = linspace(1, nHue, 300);
    plot(ax, hx, ai + ampi*cos(2*pi*(hx-1)/nHue + deg2rad(phi)), ...
        'k--', 'LineWidth', 1.5);

    peakSi = mod(-phi - shiftDeg, 360);
    txt = sprintf('peak ≈ %.0f°', peakSi);
    text(ax, 0.02, 0.95, txt, ...
        'Units','normalized','HorizontalAlignment','left', ...
        'VerticalAlignment','top');
end

% global title
if isfield(U,'exptName')
    ttl = sprintf('%s | Unit %d (%s) — Spike accounting', ...
        string(U.exptName), U.unitID, U.unitType);
else
    ttl = sprintf('%s | Unit %d (%s) — Spike accounting', ...
        string(U.dateStr), U.unitID, U.unitType);
end
sgtitle(figAcc, ttl, 'Interpreter','none');

% save
fileTag = sprintf('%s_%s_%d', string(U.dateStr), U.unitType, U.unitID);
pngName = fullfile(unitDir, sprintf('SpikeAccounting_%s.png', fileTag));
figName = fullfile(unitDir, sprintf('SpikeAccounting_%s.fig', fileTag));

exportgraphics(figAcc, pngName, 'Resolution', pngDpi);
savefig(figAcc, figName);

if ~makePlots
    close(figAcc);
end

end

% -------------------------------------------------------------------------
function [a, amp, ph, Rsq] = fit_sine_hues_local(y)
y = y(:);
n = numel(y);
h = (1:n).';
X = [ones(n,1), cos(2*pi*(h-1)/n), sin(2*pi*(h-1)/n)];
beta = X \ y;
a    = beta(1);
B    = beta(2);
C    = beta(3);
amp  = sqrt(B.^2 + C.^2);
ph   = atan2(-C, B);              % note sign to match legacy
ph   = rad2deg(ph);
yhat = X*beta;
SSres = sum((y - yhat).^2);
SStot = sum((y - mean(y)).^2);
if SStot > 0
    Rsq = 1 - SSres/SStot;
else
    Rsq = NaN;
end
end
