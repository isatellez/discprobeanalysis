function plot_boxplots(U, stats, COL, config, outDir)
% PLOT_BOXPLOTS - Boxplots of rate vs Hue per Saturation.
% FIX: Uses Synthetic Color Scaling (Vivid * Sat) for accurate RGB colors.

    if U.nTrials == 0
        return;
    end

    % --- DATA PREP ---
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
    elevID = U.trials.elevID;

    % Determine Elevation (Mode) for color lookup
    uElevs = unique(elevID(~isnan(elevID)));
    if isempty(uElevs), currentElev = 0; else, currentElev = mode(uElevs); end

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
    
    if isempty(saturIDs), maxSat = 1; else, maxSat = max(saturIDs); end

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

    % --- PLOTTING ---
    nPlotSat = nS;
    figBP = figure('Color','w','Visible',figVis,'Name','Hue box plots');
    TL = tiledlayout(figBP, 1, nPlotSat, 'Padding','compact','TileSpacing','compact');

    if isfield(U,'exptName')
        unitLabel = sprintf('%s | Unit %d (%s) — Box plots', string(U.exptName), U.unitID, U.unitType);
    else
        unitLabel = sprintf('%s | Unit %d (%s) — Box plots', string(U.dateStr), U.unitID, U.unitType);
    end
    sgtitle(figBP, unitLabel, 'Interpreter','none');

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
        
        % --- FIX: CALCULATE COLORS DYNAMICALLY ---
        % Generate color matrix for this specific saturation (nHue x 3)
        C = zeros(nHue, 3);
        for h = 1:nHue
            C(h, :) = get_synthetic_color(COL, h, sVal, maxSat, currentElev);
        end
        
        if exist('boxchart','file')
            % Modern MATLAB boxchart
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
            % Fallback for older MATLAB
            [~,~,grp] = unique(gHue);
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
    
    % per-session saves
    if ~isempty(sessionDir)
        pngName_sess = fullfile(sessionDir, sprintf('BoxPlots_%s.png', fileTag));
        figName_sess = fullfile(sessionDir, sprintf('BoxPlots_%s.fig', fileTag));
        exportgraphics(figBP, pngName_sess, 'Resolution', pngDpi);
        savefig(figBP, figName_sess);
    end
    
    % global / overall saves
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

% -------------------------------------------------------------------------
% Helper: Synthetic Color Scaling (Vivid * Saturation)
% -------------------------------------------------------------------------
function c = get_synthetic_color(COL, h, sVal, maxSat, eVal)
    % 1. Get Vivid Color (Max Saturation)
    cVivid = get_vivid_color(COL, h, maxSat, eVal);
    
    % 2. Scale brightness by Saturation factor
    if maxSat > 0
        scaleFactor = sVal / maxSat;
    else
        scaleFactor = 0;
    end
    
    c = cVivid * scaleFactor;
end

function c = get_vivid_color(COL, h, maxSat, e)
    c = [0.5 0.5 0.5]; 
    if ~isfield(COL, 'probeIDs') || ~isfield(COL, 'probeCols'), return; end
    tol = 1e-4; 
    
    % Try strict match for Max Sat
    idx = find(abs(COL.probeIDs(:,1)-h)<tol & abs(COL.probeIDs(:,2)-maxSat)<tol & abs(COL.probeIDs(:,3)-e)<tol, 1);
    
    % Fallback to Equator
    if isempty(idx)
        idx = find(abs(COL.probeIDs(:,1)-h)<tol & abs(COL.probeIDs(:,2)-maxSat)<tol & abs(COL.probeIDs(:,3)-0)<tol, 1);
    end
    
    % Fallback to Hue only (any elevation/sat if max sat missing)
    if isempty(idx)
         idx = find(abs(COL.probeIDs(:,1)-h)<tol, 1);
    end

    if ~isempty(idx)
        rawC = double(COL.probeCols(idx,:));
        if max(rawC) > 1, rawC = rawC / 255; end
        c = rawC;
    end
end