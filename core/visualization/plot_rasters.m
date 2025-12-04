function plot_rasters(U, stats, COL, config, outDir, mode)
% Wrapper around DrawRaster for DiscProbe units.
% mode: "canonical" or "mean"

if nargin < 6 || isempty(mode)
    mode = "canonical";
end
mode = string(mode);

nTr = numel(U.spk);
if nTr == 0
    return;
end

spkcell = U.spk;
mrkcell = U.mrk;

if ~iscell(spkcell) || ~iscell(mrkcell)
    error('U.spk and U.mrk should be cell arrays per trial');
end

T = U.trials;

% default grouping: all one group
groupingvar = ones(nTr,1);

% base order before we mess with it
order = (1:nTr)';

switch lower(mode)
    case "canonical"
        if all(ismember({'hueID','satID','trial'}, T.Properties.VariableNames))
            [~, order] = sortrows(T, {'hueID','satID','trial'});
        elseif all(ismember({'hueID','satID'}, T.Properties.VariableNames))
            [~, order] = sortrows(T, {'hueID','satID'});
        elseif ismember('hueID', T.Properties.VariableNames)
            [~, order] = sortrows(T, {'hueID'});
        else
            order = (1:nTr)';
        end

        if ismember('hueID', T.Properties.VariableNames)
            groupingvar = T.hueID;
            groupingvar(isnan(groupingvar)) = -1;
        else
            groupingvar = ones(nTr,1);
        end

    case "mean"
        rateEarly = compute_trial_rate(spkcell, U.winEarly);
        if isfield(U,'isOutlier') && numel(U.isOutlier) == numel(rateEarly)
            rateEarly(U.isOutlier) = NaN;
        end
        [~, order] = sort(rateEarly, 'descend', 'MissingPlacement','last');
        groupingvar = ones(nTr,1);

    otherwise
        order = (1:nTr)';
        groupingvar = ones(nTr,1);
end

% reorder everything according to the chosen order
spkcell     = spkcell(order);
mrkcell     = mrkcell(order);
groupingvar = groupingvar(order);

% per-trial colors in original trial index space, then reorder
rowcolors_all = get_trial_colors(U.trials, COL);
rowcolors     = rowcolors_all(order,:);

makePlots = true;
if isfield(config, 'plot') && isfield(config.plot, 'makePlots')
    makePlots = logical(config.plot.makePlots);
end

pngDpi = 300;
if isfield(config, 'plot') && isfield(config.plot, 'dpi')
    pngDpi = config.plot.dpi;
end

unitDir = fullfile(outDir, sprintf('%s_%d', U.unitType, U.unitID));
if ~exist(unitDir,'dir')
    mkdir(unitDir);
end

H = DrawRaster(spkcell, mrkcell, groupingvar, rowcolors);

if ~makePlots
    set(H, 'Visible', 'off');
end

ax = gca;

if isfield(U,'winFull') && numel(U.winFull) == 2
    xlim(ax, U.winFull);
end

title(ax, sprintf('%s | Unit %d (%s) — %s', ...
    string(U.dateStr), U.unitID, U.unitType, char(mode)), ...
    'Interpreter','none');

fileTag = sprintf('%s_%s_%d', string(U.dateStr), U.unitType, U.unitID);
pngName = fullfile(unitDir, sprintf('Raster_%s_%s.png', char(mode), fileTag));
figName = fullfile(unitDir, sprintf('Raster_%s_%s.fig',  char(mode), fileTag));

exportgraphics(H, pngName, 'Resolution', pngDpi);
savefig(H, figName);

if ~makePlots
    close(H);
end

end

function rate = compute_trial_rate(spkcell, win)
    dur = diff(win);
    nTr = numel(spkcell);
    rate = nan(nTr,1);
    for tt = 1:nTr
        spks = spkcell{tt};
        if isempty(spks)
            rate(tt) = 0;
        else
            spks = spks(spks >= win(1) & spks < win(2));
            rate(tt) = numel(spks) ./ dur;
        end
    end
end

function rowcolors = get_trial_colors(T, COL)
    nTr = height(T);
    rowcolors = repmat([0 0 0], nTr, 1);

    hasRGB = all(ismember({'R','G','B'}, T.Properties.VariableNames));
    if hasRGB
        rgb = double([T.R, T.G, T.B]);
        rowcolors = rgb;
        return;
    end

    if ~isfield(COL, 'probeCols') || ~isfield(COL, 'probeIDs')
        return;
    end

    ids = COL.probeIDs;
    cols = COL.probeCols;

    hueID = [];
    satID = [];
    if ismember('hueID', T.Properties.VariableNames)
        hueID = T.hueID;
    end
    if ismember('satID', T.Properties.VariableNames)
        satID = T.satID;
    end

    if isempty(hueID) || isempty(satID)
        return;
    end

    for tr = 1:nTr
        h = hueID(tr);
        s = satID(tr);

        if isnan(h) || isnan(s)
            continue;
        end

        mask = (ids(:,1) == h) & (ids(:,2) == s);
        idx = find(mask, 1, 'first');
        if ~isempty(idx)
            rowcolors(tr,:) = double(cols(idx,:));
        end
    end
end
