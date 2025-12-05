function plot_rayleigh(U, stats, COL, config, outDir)

if U.nTrials == 0
    return;
end

if ~isfield(stats, 'rayleigh') || isempty(stats.rayleigh)
    % still need trial-level info for the circle figure, so we don't bail here
end

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

% early window spikes per trial
win = U.winEarly;
dur = diff(win);

nTr = numel(U.spk);
spk_early_count = nan(nTr,1);
for tt = 1:nTr
    spks = U.spk{tt};
    if isempty(spks)
        spk_early_count(tt) = 0;
    else
        spks = spks(spks >= win(1) & spks < win(2));
        spk_early_count(tt) = numel(spks);
    end
end

if sum(spk_early_count) > 0
    w_spk = spk_early_count / sum(spk_early_count);
else
    w_spk = zeros(size(spk_early_count));
end

hueID = U.trials.hueID;
satID = U.trials.satID;

if isfield(COL, 'nHue')
    nHue = COL.nHue;
else
    nHue = max(hueID(~isnan(hueID)));
end

ZERO_HUE = 1;
if isfield(config, 'space') && isfield(config.space, 'zeroHue')
    ZERO_HUE = config.space.zeroHue;
end

thetaIdx = hueID;
thetaRad_trials = 2*pi * mod(thetaIdx - ZERO_HUE, nHue) / nHue;

hueAngles = 2*pi * mod((1:nHue) - ZERO_HUE, nHue) / nHue;
hStep = 360 / nHue;

saturIDs = unique(satID(~isnan(satID)));
saturIDs = saturIDs(:).';
nSat = numel(saturIDs);

cols_satur2 = [];
if isfield(COL, 'colsSatur')
    sz = size(COL.colsSatur);
    cols_satur2 = reshape(double(COL.colsSatur), sz(1), sz(2), 3);
end

valid = isfinite(thetaRad_trials) & isfinite(w_spk);
theta_all = thetaRad_trials(valid);
w_all     = w_spk(valid);

if isempty(theta_all) || sum(w_all)==0
    return;
end

[Rbar, mu_deg, z_all, p_all] = local_rayleigh(theta_all, w_all);

spikes_total = sum(spk_early_count);
trials_total = numel(spk_early_count);

ws_overall = normalize_to_one(sum_per_hue(thetaIdx(valid), w_all, nHue));
C_over     = get_hue_colors(cols_satur2, [], nHue);

prefHue_over = mod(round(mu_deg/hStep), nHue) + 1;

figRay = figure('Name','Rayleigh (bars)','Color','w','Visible',figVis);

if isfield(U,'exptName')
    unitTag = sprintf('%s | Unit %d (%s)', string(U.exptName), U.unitID, U.unitType);
else
    unitTag = sprintf('%s | Unit %d (%s)', string(U.dateStr), U.unitID, U.unitType);
end
sgtitle(figRay, unitTag, 'Interpreter','none');

figRay.Units = 'normalized';
AX_H   = 0.60;
AX_W   = 0.18;
BASE_Y = 0.30;
GAP_X  = 0.06;
LEFTS  = [0.05, 0.05+AX_W+GAP_X, 0.05+2*(AX_W+GAP_X), 0.05+3*(AX_W+GAP_X)];

nPlotSat = min(nSat, 3);           % only plot up to 3 saturation panels
nAxes    = 1 + nPlotSat;           % overall + per-sat

pax = gobjects(1, nAxes);
for i = 1:nAxes
    pax(i) = polaraxes('Parent',figRay, ...
                       'Units','normalized', ...
                       'Position',[LEFTS(i) BASE_Y AX_W AX_H]);
end

SPOKE_LW = 3.5;
MU_LW    = 3.0;
STATS_FS = 11;
ANG_EPS  = 1e-3;

ax0 = pax(1);
title(ax0, 'Overall', 'FontWeight','bold');
hold(ax0,'on');
style_polar(ax0, nHue, false);

draw_spoke_bars_line(ax0, ws_overall, hueAngles, C_over, SPOKE_LW);

mu = deg2rad(mu_deg);
polarplot(ax0, [mu-ANG_EPS mu+ANG_EPS], [0 Rbar], '-', 'LineWidth', MU_LW, 'Color','k');
polarplot(ax0, mu, Rbar, 'ko', 'MarkerFaceColor','k', 'MarkerSize',5);

txt_over = sprintf('R=%.3f, p=%.5g, z=%.2f, \\mu=%0.0f^\\circ (hue %d)\nspikes=%d, trials=%d', ...
    Rbar, p_all, z_all, mu_deg, prefHue_over, spikes_total, trials_total);
place_stats_below(figRay, ax0, txt_over, STATS_FS);

R_s  = nan(1,nPlotSat);
mu_s = nan(1,nPlotSat);
z_s  = nan(1,nPlotSat);
p_s  = nan(1,nPlotSat);
tr_s = nan(1,nPlotSat);
sp_s = nan(1,nPlotSat);

for si = 1:nPlotSat
    ax = pax(si+1);
    sVal = saturIDs(si);

    title(ax, sprintf('sat = %.2f', sVal), 'FontWeight','bold');
    hold(ax,'on');
    style_polar(ax, nHue, false);

    sel  = satID == sVal & isfinite(thetaRad_trials) & isfinite(w_spk);
    th_i = thetaRad_trials(sel);
    w_i  = w_spk(sel);

    tr_s(si) = sum(sel);
    sp_s(si) = sum(spk_early_count(sel));

    if isempty(th_i) || sum(w_i)==0
        Rsi   = 0;
        mu_si = 0;
        z_si  = NaN;
        p_si  = NaN;
        ws_si = zeros(nHue,1);
    else
        [Rsi, mu_si, z_si, p_si] = local_rayleigh(th_i, w_i);
        ws_si = normalize_to_one(sum_per_hue(thetaIdx(sel), w_i, nHue));
    end

    R_s(si)  = Rsi;
    mu_s(si) = mu_si;
    z_s(si)  = z_si;
    p_s(si)  = p_si;

    C_si = get_hue_colors(cols_satur2, si, nHue);
    draw_spoke_bars_line(ax, ws_si, hueAngles, C_si, SPOKE_LW);

    mu = deg2rad(mu_si);
    polarplot(ax,[mu-ANG_EPS mu+ANG_EPS],[0 Rsi],'-','LineWidth',MU_LW,'Color','k');
    polarplot(ax,mu,Rsi,'ko','MarkerFaceColor','k','MarkerSize',5);

    prefHue_si = mod(round(mu_si/hStep), nHue)+1;
    txt_si = sprintf('R=%.3f, p=%.5g, z=%.2f, \\mu=%0.0f^\\circ (hue %d)\nspikes=%d, trials=%d', ...
        Rsi, p_si, z_si, mu_si, prefHue_si, sp_s(si), tr_s(si));
    place_stats_below(figRay, ax, txt_si, STATS_FS);
end

fileTag = sprintf('%s_%s_%d', string(U.dateStr), U.unitType, U.unitID);
pngName = fullfile(unitDir, sprintf('Rayleigh_%s.png', fileTag));
figName = fullfile(unitDir, sprintf('Rayleigh_%s.fig', fileTag));

exportgraphics(figRay, pngName, 'Resolution', pngDpi);
savefig(figRay, figName);

if ~makePlots
    close(figRay);
end

end

% helpers

function [Rbar, mu_deg, z, p] = local_rayleigh(thetaRad, w)
thetaRad = thetaRad(:);
w        = w(:);
C = sum(w .* cos(thetaRad));
S = sum(w .* sin(thetaRad));
W = sum(w);
Rbar = hypot(C,S) / max(W, eps);
mu   = atan2(S, C);
mu_deg = mod(rad2deg(mu), 360);
z = W * Rbar.^2;
p = exp(-z);
end

function x = normalize_to_one(x)
x = x(:);
if all(~isfinite(x)) || all(x==0 | isnan(x))
    x(:) = 0;
    return;
end
x(x < 0) = 0;
m = max(x);
if m > 0
    x = x ./ m;
end
end

function ws = sum_per_hue(hueIdx, w, nHue)
hueIdx = hueIdx(:);
w      = w(:);
valid  = isfinite(hueIdx) & isfinite(w) & hueIdx>=1 & hueIdx<=nHue;
ws = accumarray(hueIdx(valid), w(valid), [nHue 1], @sum, 0);
end

function C = get_hue_colors(cols_satur2, satIdx, nHue)
if isempty(cols_satur2)
    C = repmat([0.5 0.5 0.5], nHue, 1);
    return;
end

sz = size(cols_satur2);
if numel(sz) ~= 3
    C = repmat([0.5 0.5 0.5], nHue, 1);
    return;
end

if isempty(satIdx) || isnan(satIdx)
    satIdx = sz(1);
end
satIdx = min(max(satIdx,1), sz(1));

C = squeeze(cols_satur2(satIdx, 1:nHue, :));
if size(C,1) ~= nHue
    C = C(1:min(end,nHue),:);
end
C = double(C) ./ 255;
end

function style_polar(ax, nHue, showRadialLabels) %#ok<INUSD>
if nargin < 3
    showRadialLabels = false;
end

thetaticks(ax, 0:45:315);
rticks(ax, 0:0.5:1);
rlim(ax,[0 1]);
if ~showRadialLabels
    ax.RTickLabel = {};
end
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
end

function draw_spoke_bars_line(ax, ws, hueAngles, C, lw)
ws = ws(:);
hueAngles = hueAngles(:);
n = numel(ws);
for k = 1:n
    r = ws(k);
    if r <= 0, continue; end
    ang = hueAngles(k);
    c = C(min(k, size(C,1)),:);
    polarplot(ax, [ang ang], [0 r], '-', 'LineWidth', lw, 'Color', c);
end
end

function place_stats_below(figH, ax, txt, fs)
axPos = ax.Position;
x = axPos(1) + axPos(3)/2;
y = axPos(2) - 0.06;

annotation(figH, 'textbox', [x-0.12 y 0.24 0.08], ...
    'String', txt, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','top', ...
    'EdgeColor','none', ...
    'FontSize', fs, ...
    'Interpreter','tex');
end
