function plot_boxplots(U, stats, COL, config, outDir)

if U.nTrials == 0
    return;
end

% early window, same as detect_outliers
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

hueID = U.trials.hueID;
satID = U.trials.satID;

% defaults matching original
USE_NORM = false;
METHOD   = 'tukey';
MIN_N    = 5;

if isfield(config,'boxPlots')
    bp = config.boxPlots;
    if isfield(bp,'useNorm'), USE_NORM = logical(bp.useNorm); end
    if isfield(bp,'method'),  METHOD   = lower(bp.method);    end
    if isfield(bp,'minN'),    MIN_N    = bp.minN;             end
end

if isfield(COL,'nHue')
    nHue = COL.nHue;
else
    nHue = max(hueID(~isnan(hueID)));
end

ZERO_HUE = 1;
if isfield(config,'space') && isfield(config.space,'zeroHue')
    ZERO_HUE = config.space.zeroHue;
end

% saturation list and indices per trial
saturIDs = unique(satID(~isnan(satID)));
saturIDs = saturIDs(:).';
nS       = numel(saturIDs);

[~, satIdx_u] = ismember(satID, saturIDs);
hueIdx_u = hueID;

% choose values and y label
if USE_NORM
    r = rateEarly;
    r = r - min(r);
    denom = max(r);
    if denom <= 0
        vals = zeros(size(r));
    else
        vals = r ./ denom;
    end
    ylab = 'Normalized rate (0–1)';
else
    vals = rateEarly;
    ylab = 'Rate (Hz)';
end

% figure visibility toggle
makePlots = true;
if isfield(config,'plot') && isfield(config.plot,'makePlots')
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

Upaths  = get_unit_paths(config, dateStr, unitType, unitID);
unitDir = Upaths.figures;
if ~exist(unitDir,'dir')
    mkdir(unitDir);
end


% per-session boxplots dir under figs/discprobe/boxplots
if nargin < 5 || isempty(outDir)
    outDir = '';
end
sessionDir = '';
if ~isempty(outDir)
    sessionDir = fullfile(outDir, 'boxplots');
    if ~exist(sessionDir,'dir')
        mkdir(sessionDir);
    end
end

% global "overall" dir under output/figs/discprobe/boxplots
globalDir = '';
if isfield(config,'paths') && isfield(config.paths,'globalDiscProbeFigRoot') ...
        && ~isempty(config.paths.globalDiscProbeFigRoot)
    globalDir = fullfile(config.paths.globalDiscProbeFigRoot, 'boxplots');
    if ~exist(globalDir,'dir')
        mkdir(globalDir);
    end
end

% one subplot per saturation
nPlotSat = nS;
figBP = figure('Color','w','Visible',figVis,'Name','Hue box plots');
TL = tiledlayout(figBP, 1, nPlotSat, 'Padding','compact','TileSpacing','compact');

if isfield(U,'exptName')
    unitLabel = sprintf('%s | Unit %d (%s) — Box plots', string(U.exptName), U.unitID, U.unitType);
else
    unitLabel = sprintf('%s | Unit %d (%s) — Box plots', string(U.dateStr), U.unitID, U.unitType);
end
sgtitle(figBP, unitLabel, 'Interpreter','none');

% color table from COL.colsSatur
cols_satur2 = [];
if isfield(COL,'colsSatur') && ~isempty(COL.colsSatur)
    sz = size(COL.colsSatur);                         % [nSat x nHue x 3]
    cols_satur2 = reshape(double(COL.colsSatur), sz(1), sz(2), 3);
end

for si = 1:nPlotSat
    ax = nexttile(TL);
    hold(ax,'on');

    sVal = saturIDs(si);
    msk  = (satIdx_u == si) & isfinite(hueIdx_u) & isfinite(vals);

    if ~any(msk)
        title(ax, sprintf('sat=%.2f (no data)', sVal));
        axis(ax,'off');
        continue;
    end

    gHue = hueIdx_u(msk);
    x    = vals(msk);

    if exist('boxchart','file')
        % build hue-wise colors for this sat
        if ~isempty(cols_satur2)
            nSatCols = size(cols_satur2,1);
            satRow   = min(si, nSatCols);
            C = squeeze(cols_satur2(satRow,1:nHue,:));   % nHue×3
            C = double(C);
            if max(C(:)) > 1, C = C/255; end
            C(~isfinite(C)) = 0;
        else
            C = repmat([0.5 0.5 0.5], nHue, 1);
        end

        for h = 1:nHue
            j = find(msk & hueIdx_u == h);
            if isempty(j), continue; end

            b = boxchart(hueIdx_u(j), vals(j), 'Parent', ax, ...
                'MarkerStyle','.', 'BoxFaceAlpha',0.6);

            rgb = C(h,:);
            rgb = min(max(rgb,0),1);
            b.BoxFaceColor = rgb;
            b.BoxEdgeColor = rgb;
            b.MarkerColor  = rgb;
        end
    else
        % fallback to boxplot with per-hue grouping
        g = gHue;
        C = repmat([0 0 0], nHue, 1);
        if ~isempty(cols_satur2)
            nSatCols = size(cols_satur2,1);
            satRow   = min(si, nSatCols);
            C = squeeze(cols_satur2(satRow,1:nHue,:));
            C = double(C);
            if max(C(:)) > 1, C = C/255; end
            C(~isfinite(C)) = 0;
        end
        [~,~,grp] = unique(g);
        boxplot(ax, x, grp, 'Symbol','.', 'Colors', C(1:max(grp),:));
        set(ax,'XTick',1:nHue,'XLim',[0.5 nHue+0.5]);
    end

    title(ax, sprintf('sat=%.2f', sVal));
    xlabel(ax,'Hue');
    ylabel(ax, ylab);

    % label x-axis as degrees if DKL
    if isfield(config,'space') && isfield(config.space,'mode') && strcmpi(config.space.mode,'dkl')
        degLabels = mod((1:nHue) - ZERO_HUE, nHue) * (360/nHue);
        set(ax,'XTick',1:nHue,'XTickLabel', ...
            arrayfun(@(d)sprintf('%d',round(d)),degLabels,'UniformOutput',false));
        xlabel(ax,'Angle (deg)');
    else
        set(ax,'XTick',1:nHue);
    end

    ax.XLim = [0.5, nHue+0.5];
    if USE_NORM, ylim(ax,[0 1]); end
    grid(ax,'on');

    % outliers overlay
    if isfield(U,'isOutlier') && numel(U.isOutlier) == nTr
        out_msk = U.isOutlier & msk;
        if any(out_msk)
            plot(ax, hueIdx_u(out_msk), vals(out_msk), 'ro', 'MarkerSize',5, 'LineWidth',1.2);
        end
    end
end

pngDpi = 200;
if isfield(config,'plot') && isfield(config.plot,'dpi')
    pngDpi = config.plot.dpi;
end

fileTag = sprintf('%s_%s_%d', string(U.dateStr), U.unitType, U.unitID);

% per-unit saves
pngName_unit = fullfile(unitDir, sprintf('BoxPlots_%s.png', fileTag));
figName_unit = fullfile(unitDir, sprintf('BoxPlots_%s.fig', fileTag));
exportgraphics(figBP, pngName_unit, 'Resolution', pngDpi);
savefig(figBP, figName_unit);

% per-session saves (if sessionDir defined)
if ~isempty(sessionDir)
    pngName_sess = fullfile(sessionDir, sprintf('BoxPlots_%s.png', fileTag));
    figName_sess = fullfile(sessionDir, sprintf('BoxPlots_%s.fig', fileTag));
    exportgraphics(figBP, pngName_sess, 'Resolution', pngDpi);
    savefig(figBP, figName_sess);
end

% global / overall saves (if globalDir defined)
if ~isempty(globalDir)
    pngName_glob = fullfile(globalDir, sprintf('BoxPlots_%s.png', fileTag));
    figName_glob = fullfile(globalDir, sprintf('BoxPlots_%s.fig', fileTag));
    exportgraphics(figBP, pngName_glob, 'Resolution', pngDpi);
    savefig(figBP, figName_glob);
end

if ~makePlots
    close(figBP);
end

end
