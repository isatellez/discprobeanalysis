function plot_rayleigh(U, stats, COL, config, outDir)
% PLOT_RAYLEIGH - Polar tuning plots.
% UPDATES:
% 1. Orientation: 0° Right, 90° Top (Counter-Clockwise).
% 2. Alignment: Angles calculated relative to config.space.zeroHue.
% 3. Cleanliness: Axis lines and grids are hidden (only data shown).

    if U.nTrials == 0
        return;
    end
    
    % Check stats
    if ~isfield(stats, 'rayleigh') || isempty(stats.rayleigh)
        % proceed anyway
    end
    
    % --- CONFIG ---
    makePlots = true;
    if isfield(config, 'plot') && isfield(config.plot, 'makePlots')
        makePlots = logical(config.plot.makePlots);
    end
    figVis = 'on';
    if ~makePlots, figVis = 'off'; end
    
    pngDpi = 300;
    if isfield(config, 'plot') && isfield(config.plot, 'dpi')
        pngDpi = config.plot.dpi;
    end
    
    doSaveFig = true;
    if isfield(config, 'plot') && isfield(config.plot, 'saveFigs')
        doSaveFig = config.plot.saveFigs;
    end
    
    % --- PATHS ---
    dateStr  = char(string(U.dateStr));
    unitType = char(string(U.unitType));
    unitID   = U.unitID;
    
    Upaths  = get_unit_paths(config, dateStr, unitType, unitID);
    unitDir = Upaths.figures;
    if ~exist(unitDir,'dir'), mkdir(unitDir); end
    
    if nargin < 5 || isempty(outDir), outDir = ''; end
    sessionDir = '';
    if ~isempty(outDir)
        useSub = true;
        if isfield(config, 'plot') && isfield(config.plot, 'suppressSubfolders') && config.plot.suppressSubfolders
            useSub = false;
        end
        if useSub
            sessionDir = fullfile(outDir, 'rayleigh');
        else
            sessionDir = outDir; 
        end
        if ~exist(sessionDir,'dir'), mkdir(sessionDir); end
    end
    
    globalDir = '';
    if isfield(config,'paths') && isfield(config.paths,'globalDiscProbeFigRoot') ...
            && ~isempty(config.paths.globalDiscProbeFigRoot)
        globalDir = fullfile(config.paths.globalDiscProbeFigRoot, 'rayleigh');
        if ~exist(globalDir,'dir'), mkdir(globalDir); end
    end
    
    % --- DATA PREP ---
    win = U.winEarly;
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
    
    hueID  = U.trials.hueID;
    satID  = U.trials.satID;
    elevID = U.trials.elevID;
    
    uElevs = unique(elevID(~isnan(elevID)));
    if isempty(uElevs), currentElev = 0; else, currentElev = mode(uElevs); end
    
    if isfield(COL, 'nHue')
        nHue = COL.nHue;
    else
        nHue = max(hueID(~isnan(hueID)));
    end
    
    % --- ANGLE ALIGNMENT ---
    % Default to 16 if not specified, per your request
    ZERO_HUE = 16; 
    if isfield(config, 'space') && isfield(config.space, 'zeroHue')
        ZERO_HUE = config.space.zeroHue;
    end
    
    thetaIdx        = hueID;
    % Calculate radians relative to Zero Hue
    % e.g. if Zero=16, then Hue 16 -> 0 rad. Hue 4 -> pi/2 rad (90 deg)
    thetaRad_trials = 2*pi * mod(thetaIdx - ZERO_HUE, nHue) / nHue;
    
    % Pre-calculate angles for drawing bars
    hueAngles       = 2*pi * mod((1:nHue) - ZERO_HUE, nHue) / nHue;
    hStep           = 360 / nHue;
    
    saturIDs = unique(satID(~isnan(satID)));
    saturIDs = saturIDs(:).';
    nSat     = numel(saturIDs);
    if isempty(saturIDs), maxSat = 1; else, maxSat = max(saturIDs); end
    
    valid     = isfinite(thetaRad_trials) & isfinite(w_spk);
    theta_all = thetaRad_trials(valid);
    w_all     = w_spk(valid);
    
    if isempty(theta_all) || sum(w_all) == 0, return; end
    
    % Statistics
    [Rbar, mu_deg, z_all, p_all] = local_rayleigh(theta_all, w_all);
    
    spikes_total = sum(spk_early_count);
    trials_total = numel(spk_early_count);
    
    ws_overall   = normalize_to_one(sum_per_hue(thetaIdx(valid), w_all, nHue));
    prefHue_over = mod(round(mu_deg/hStep), nHue) + 1;
    
    % --- Overall Color Map (Full Vivid) ---
    C_over = zeros(nHue, 3);
    for h = 1:nHue
        C_over(h, :) = get_synthetic_color(COL, h, maxSat, maxSat, currentElev);
    end
    
    % ---------- Figure Layout ----------
    figRay = figure('Name','Rayleigh (bars)','Color','w','Visible',figVis);
    
    nCols = min(nSat, 4); if nCols < 1, nCols = 1; end
    nRowsSat = ceil(nSat / nCols);
    totalRows = 1 + nRowsSat;
    
    t = tiledlayout(figRay, totalRows, nCols, 'TileSpacing','normal', 'Padding','compact');

    if isfield(U,'exptName')
        unitTag = sprintf('%s | Unit %d (%s)', string(U.exptName), U.unitID, U.unitType);
    else
        unitTag = sprintf('%s | Unit %d (%s)', string(U.dateStr), U.unitID, U.unitType);
    end
    sgtitle(figRay, unitTag, 'Interpreter','none');
    
    figHeight = 350 + (nRowsSat * 300); 
    figWidth  = 300 * nCols;
    figRay.Position = [100, 100, min(figWidth, 1600), min(figHeight, 1000)];

    SPOKE_LW = 3.5;
    MU_LW    = 3.0;
    ANG_EPS  = 1e-3;
    
    % --- Plot 1: Overall ---
    ax0 = nexttile(t, 1, [1, nCols]); 
    ax0 = polaraxes('Parent', t);
    ax0.Layout.Tile = 1;
    ax0.Layout.TileSpan = [1, nCols];
    
    title(ax0, 'Overall', 'FontWeight','bold');
    hold(ax0,'on');
    style_polar(ax0);
    
    draw_spoke_bars_line(ax0, ws_overall, hueAngles, C_over, SPOKE_LW);
    
    mu = deg2rad(mu_deg);
    polarplot(ax0, [mu-ANG_EPS mu+ANG_EPS], [0 Rbar], '-', 'LineWidth', MU_LW, 'Color','k');
    polarplot(ax0, mu, Rbar, 'ko', 'MarkerFaceColor','k', 'MarkerSize',5);
    
    txt_over = sprintf('R=%.3f, p=%.4g, z=%.1f, \\mu=%0.0f^\\circ (hue %d)\nspikes=%d, trials=%d', ...
        Rbar, p_all, z_all, mu_deg, prefHue_over, spikes_total, trials_total);
    title(ax0, { 'Overall', txt_over }, 'FontWeight','bold', 'FontSize', 9);
    
    % --- Per-Sat Panels ---
    for si = 1:nSat
        ax = nexttile(t);
        ax = polaraxes('Parent', t);
        ax.Layout.Tile = nCols + si; 
        
        sVal = saturIDs(si);
        hold(ax,'on');
        style_polar(ax);
    
        sel  = satID == sVal & isfinite(thetaRad_trials) & isfinite(w_spk);
        th_i = thetaRad_trials(sel);
        w_i  = w_spk(sel);
        
        tr_si_count = sum(sel);
        sp_si_count = sum(spk_early_count(sel));
    
        if isempty(th_i) || sum(w_i)==0
            Rsi = 0; mu_si = 0; z_si = NaN; p_si = NaN; ws_si = zeros(nHue,1);
        else
            [Rsi, mu_si, z_si, p_si] = local_rayleigh(th_i, w_i);
            ws_si = normalize_to_one(sum_per_hue(thetaIdx(sel), w_i, nHue));
        end
    
        % --- Synthetic Color for this Saturation ---
        C_si = zeros(nHue, 3);
        for h = 1:nHue
            C_si(h, :) = get_synthetic_color(COL, h, sVal, maxSat, currentElev);
        end
        
        draw_spoke_bars_line(ax, ws_si, hueAngles, C_si, SPOKE_LW);
        
        mu = deg2rad(mu_si);
        polarplot(ax,[mu-ANG_EPS mu+ANG_EPS],[0 Rsi],'-','LineWidth',MU_LW,'Color','k');
        polarplot(ax,mu,Rsi,'ko','MarkerFaceColor','k','MarkerSize',5);
        
        txt_si = sprintf('R=%.2f, p=%.3g\n\\mu=%0.0f^\\circ, n=%d', Rsi, p_si, mu_si, sp_si_count);
        title(ax, { sprintf('sat = %.2f', sVal), txt_si }, 'FontWeight','normal', 'FontSize', 8);
    end
    
    % ---------- Saving ----------
    fileTag  = sprintf('%s_%s_%d', string(U.dateStr), U.unitType, U.unitID);
    baseName = sprintf('Rayleigh_%s', fileTag);
    
    exportgraphics(figRay, fullfile(unitDir, [baseName '.png']), 'Resolution', pngDpi);
    if doSaveFig, savefig(figRay, fullfile(unitDir, [baseName '.fig'])); end
    
    if ~isempty(sessionDir)
        exportgraphics(figRay, fullfile(sessionDir, [baseName '.png']), 'Resolution', pngDpi);
        if doSaveFig, savefig(figRay, fullfile(sessionDir, [baseName '.fig'])); end
    end
    
    if ~isempty(globalDir)
        exportgraphics(figRay, fullfile(globalDir, [baseName '.png']), 'Resolution', pngDpi);
        if doSaveFig, savefig(figRay, fullfile(globalDir, [baseName '.fig'])); end
    end
    
    if ~makePlots, close(figRay); end
end

% ---------- Helpers ----------

function [Rbar, mu_deg, z, p] = local_rayleigh(thetaRad, w)
    thetaRad = thetaRad(:); w = w(:);
    C = sum(w .* cos(thetaRad)); S = sum(w .* sin(thetaRad)); W = sum(w);
    Rbar = hypot(C,S) / max(W, eps); mu = atan2(S, C); mu_deg = mod(rad2deg(mu), 360);
    z = W * Rbar.^2; p = exp(-z);
end

function x = normalize_to_one(x)
    x = x(:); if all(~isfinite(x)) || all(x==0 | isnan(x)), x(:)=0; return; end
    x(x < 0) = 0; m = max(x); if m > 0, x = x ./ m; end
end

function ws = sum_per_hue(hueIdx, w, nHue)
    hueIdx = hueIdx(:); w = w(:); valid = isfinite(hueIdx) & isfinite(w) & hueIdx>=1 & hueIdx<=nHue;
    ws = accumarray(hueIdx(valid), w(valid), [nHue 1], @sum, 0);
end

% -------------------------------------------------------------------------
% POLAR STYLE: DKL Space (0 Right, 90 Top, No Axes)
% -------------------------------------------------------------------------
function style_polar(ax)
    % Orientation
    ax.ThetaZeroLocation = 'right';
    ax.ThetaDir = 'counterclockwise';
    
    % Remove axis lines, grids, and ticks for clean look
    ax.ThetaAxis.Visible = 'off';
    ax.RAxis.Visible = 'off';
    ax.GridColor = 'none';
    ax.GridAlpha = 0;
    
    rlim(ax, [0 1]);
end

function draw_spoke_bars_line(ax, ws, hueAngles, C, lw)
    ws = ws(:); hueAngles = hueAngles(:);
    for k = 1:numel(ws)
        if ws(k) <= 0, continue; end
        c = C(min(k, size(C,1)),:);
        polarplot(ax, [hueAngles(k) hueAngles(k)], [0 ws(k)], '-', 'LineWidth', lw, 'Color', c);
    end
end

% -------------------------------------------------------------------------
% Helper: Synthetic Color Scaling (Guarantees consistency)
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
    idx = find(abs(COL.probeIDs(:,1)-h)<tol & abs(COL.probeIDs(:,2)-maxSat)<tol & abs(COL.probeIDs(:,3)-e)<tol, 1);
    if isempty(idx), idx = find(abs(COL.probeIDs(:,1)-h)<tol & abs(COL.probeIDs(:,2)-maxSat)<tol & abs(COL.probeIDs(:,3)-0)<tol, 1); end
    if isempty(idx), idx = find(abs(COL.probeIDs(:,1)-h)<tol, 1); end
    if ~isempty(idx), rawC = double(COL.probeCols(idx,:)); if max(rawC)>1, rawC=rawC/255; end; c=rawC; end
end