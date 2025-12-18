function plot_rasters(U, stats, COL, config, outDir, mode)

if U.nTrials == 0
    return;
end

if nargin < 6 || isempty(mode)
    mode = "canonical";
end
mode = string(mode);

% figure visibility
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

% ---------- paths: per-unit, per-session, global ----------

dateStr  = char(string(U.dateStr));
unitType = char(string(U.unitType));
unitID   = U.unitID;

% per-unit dir: units/<unit>/figures (flat, no subfolders)
Upaths  = get_unit_paths(config, dateStr, unitType, unitID);
unitDir = Upaths.figures;
if ~exist(unitDir,'dir')
    mkdir(unitDir);
end


% per-session dir: <date>/figs/discprobe/rasters
if nargin < 5 || isempty(outDir)
    outDir = '';
end
sessionDir = '';
if ~isempty(outDir)
    sessionDir = fullfile(outDir, 'rasters');
    if ~exist(sessionDir,'dir')
        mkdir(sessionDir);
    end
end

% global dir: output/figs/discprobe/rasters
globalDir = '';
if isfield(config,'paths') && isfield(config.paths,'globalDiscProbeFigRoot') ...
        && ~isempty(config.paths.globalDiscProbeFigRoot)
    globalDir = fullfile(config.paths.globalDiscProbeFigRoot, 'rasters');
    if ~exist(globalDir,'dir')
        mkdir(globalDir);
    end
end

% ---------- data prep ----------

spkcell_u = U.spk(:);
mrkcell_u = U.mrk(:);

if isempty(spkcell_u)
    return;
end

% row colors per trial (pastel bands)
rowRGB = [];
if isfield(U.trials,'rowRGB')
    rowRGB = double(U.trials.rowRGB);
elseif all(isfield(U.trials,{'R','G','B'}))
    rowRGB = double([U.trials.R, U.trials.G, U.trials.B]);
end
if ~isempty(rowRGB)
    if max(rowRGB,[],'all') > 1
        rowRGB = rowRGB ./ 255;
    end
else
    if isfield(COL,'colsRGB') && ~isempty(COL.colsRGB)
        cols = double(COL.colsRGB);
        if max(cols(:))>1, cols = cols/255; end
        nT   = numel(spkcell_u);
        rowRGB = cols(mod((1:nT)-1,size(cols,1))+1, :);
    else
        rowRGB = repmat([0.9 0.9 0.9], numel(spkcell_u), 1);
    end
end

% trial-type grouping variable
if isfield(U.trials,'trialTypeID')
    gid_base = U.trials.trialTypeID(:);
elseif isfield(U.trials,'trType')
    gid_base = U.trials.trType(:,1);
elseif isfield(U.trials,'hueID')
    gid_base = U.trials.hueID(:);
else
    gid_base = (1:U.nTrials).';
end

% early window rate per trial
win = U.winEarly;
dur = diff(win);
nTr = numel(spkcell_u);
spk_early_hz = nan(nTr,1);
for tt = 1:nTr
    sp = spkcell_u{tt};
    if isempty(sp)
        spk_early_hz(tt) = 0;
    else
        sp = sp(sp >= win(1) & sp < win(2));
        spk_early_hz(tt) = numel(sp) / dur;
    end
end

nTrials_now = numel(spkcell_u);
nNonzero    = sum(cellfun(@(x) ~isempty(x), spkcell_u));
nSpikes_tot = sum(cellfun(@numel, spkcell_u));

% labels
if strcmpi(U.unitType,'SU')
    unitLabel = sprintf('UnitID %d (SU)', U.unitID);
else
    unitLabel = sprintf('MUA ID %d', U.unitID);
end

if isfield(U,'unitIdx')
    unitNumStr = sprintf('UnitNUM %d', U.unitIdx);
else
    unitNumStr = '';
end

sessionStamp = string(U.dateStr);
fileBaseTag  = sprintf('%s_%s_%d', string(U.dateStr), U.unitType, U.unitID);

% ---------- plotting modes ----------

switch mode
    case "canonical"
        groupvar = gid_base;

        H = DrawRaster(spkcell_u, mrkcell_u, groupvar, rowRGB);
        set(H,'Color','w','Visible',figVis,'Name','Raster — canonical');

        if ~isempty(unitNumStr)
            ttl = sprintf('%s | %s — %s | Sorted by trial type order | trials=%d (nonzero=%d) | spikes=%d', ...
                sessionStamp, unitLabel, unitNumStr, nTrials_now, nNonzero, nSpikes_tot);
        else
            ttl = sprintf('%s | %s | Sorted by trial type order | trials=%d (nonzero=%d) | spikes=%d', ...
                sessionStamp, unitLabel, nTrials_now, nNonzero, nSpikes_tot);
        end
        title(ttl,'Interpreter','none','FontWeight','bold');

        baseName = sprintf('Raster_canonical_%s', fileBaseTag);

        % per-unit
        png_unit = fullfile(unitDir, [baseName '.png']);
        fig_unit = fullfile(unitDir, [baseName '.fig']);
        exportgraphics(H, png_unit, 'Resolution', pngDpi);
        savefig(H, fig_unit);

        % per-session
        if ~isempty(sessionDir)
            png_sess = fullfile(sessionDir, [baseName '.png']);
            fig_sess = fullfile(sessionDir, [baseName '.fig']);
            exportgraphics(H, png_sess, 'Resolution', pngDpi);
            savefig(H, fig_sess);
        end

        % global
        if ~isempty(globalDir)
            png_glob = fullfile(globalDir, [baseName '.png']);
            fig_glob = fullfile(globalDir, [baseName '.fig']);
            exportgraphics(H, png_glob, 'Resolution', pngDpi);
            savefig(H, fig_glob);
        end

        if ~makePlots && ishghandle(H)
            close(H);
        end

    case "mean"
        % group by trial type, sort by mean early rate
        gid = gid_base;
        [ug,~,gix] = unique(gid,'stable');
        m_by_group = accumarray(gix, spk_early_hz, [numel(ug) 1], @mean, NaN);

        [~, gorder] = sort(m_by_group, 'descend', 'MissingPlacement','last');

        ix = cell2mat(arrayfun(@(k) find(gix==gorder(k)).', 1:numel(gorder), 'UniformOutput', false));

        rank_of_group         = zeros(numel(ug),1);
        rank_of_group(gorder) = 1:numel(gorder);
        grp_sort              = rank_of_group(gix(ix));

        spk_sort = spkcell_u(ix);
        mrk_sort = mrkcell_u(ix);
        rgb_sort = rowRGB(ix,:);

        H2 = DrawRaster(spk_sort, mrk_sort, grp_sort, rgb_sort);
        set(H2,'Color','w','Visible',figVis,'Name','Raster — mean spikes');

        if ~isempty(unitNumStr)
            ttl = sprintf('%s | %s — %s | Sorted by mean spikes (high→low) | trials=%d (nonzero=%d) | spikes=%d', ...
                sessionStamp, unitLabel, unitNumStr, nTrials_now, nNonzero, nSpikes_tot);
        else
            ttl = sprintf('%s | %s | Sorted by mean spikes (high→low) | trials=%d (nonzero=%d) | spikes=%d', ...
                sessionStamp, unitLabel, nTrials_now, nNonzero, nSpikes_tot);
        end
        title(ttl,'Interpreter','none','FontWeight','bold');

        baseName = sprintf('Raster_byGroupMean_%s', fileBaseTag);

        % per-unit
        png_unit = fullfile(unitDir, [baseName '.png']);
        fig_unit = fullfile(unitDir, [baseName '.fig']);
        exportgraphics(H2, png_unit, 'Resolution', pngDpi);
        savefig(H2, fig_unit);

        % per-session
        if ~isempty(sessionDir)
            png_sess = fullfile(sessionDir, [baseName '.png']);
            fig_sess = fullfile(sessionDir, [baseName '.fig']);
            exportgraphics(H2, png_sess, 'Resolution', pngDpi);
            savefig(H2, fig_sess);
        end

        % global
        if ~isempty(globalDir)
            png_glob = fullfile(globalDir, [baseName '.png']);
            fig_glob = fullfile(globalDir, [baseName '.fig']);
            exportgraphics(H2, png_glob, 'Resolution', pngDpi);
            savefig(H2, fig_glob);
        end

        if ~makePlots && ishghandle(H2)
            close(H2);
        end

    otherwise
        warning('plot_rasters: unknown mode "%s"', mode);
end

end
