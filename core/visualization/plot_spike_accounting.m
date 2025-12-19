function plot_spike_accounting(U, stats, COL, config, outDir, varargin)
% PLOT_SPIKE_ACCOUNTING
% FIXES: 
% 1. Stats: Prints R^2, Amp, Peak on all panels.
% 2. Labels: X-axis shows sparse DKL angles (0, 90, etc) instead of Hue Index.
% 3. Colors: Uses Synthetic Scaling (VividColor * Saturation).

    if U.nTrials == 0
        return;
    end
    
    % --- 1. CONFIG CHECK ---
    makePlots = true;
    if isfield(config,'plot') && isfield(config.plot,'makePlots')
        makePlots = logical(config.plot.makePlots);
    end
    
    figVis = 'on';
    if ~makePlots, figVis = 'off'; end
    
    pngDpi = 300;
    if isfield(config,'plot') && isfield(config.plot,'dpi')
        pngDpi = config.plot.dpi;
    end
    
    doSaveFig = true;
    if isfield(config, 'plot') && isfield(config.plot, 'saveFigs')
        doSaveFig = config.plot.saveFigs;
    end

    sessionSubdir = 'spikeaccounting';
    if nargin >= 6 && ~isempty(varargin{1})
        sessionSubdir = varargin{1};
    end
    
    % --- 2. PATH SETUP ---
    dateStr  = char(string(U.dateStr));
    unitType = char(string(U.unitType));
    unitID   = U.unitID;
    
    Upaths  = get_unit_paths(config, dateStr, unitType, unitID);
    unitDir = Upaths.figures;
    if ~exist(unitDir,'dir'), mkdir(unitDir); end
    
    if nargin < 5 || isempty(outDir), outDir = ''; end
    sessionDir = '';
    if ~isempty(outDir)
        if contains(outDir, 'sphere_slices')
            sessionDir = outDir;
        else
            if isempty(sessionSubdir)
                sessionDir = outDir;
            else
                sessionDir = fullfile(outDir, sessionSubdir);
            end
        end
        if ~exist(sessionDir,'dir'), mkdir(sessionDir); end
    end
    
    globalDir = '';
    if isfield(config,'paths') && isfield(config.paths,'globalDiscProbeFigRoot') ...
            && ~isempty(config.paths.globalDiscProbeFigRoot)
        globalDir = fullfile(config.paths.globalDiscProbeFigRoot, 'spikeaccounting');
        if ~exist(globalDir,'dir'), mkdir(globalDir); end
    end
    
    % --- 3. DATA PROCESSING ---
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
    elevID = U.trials.elevID;
    
    uElevs = unique(elevID(~isnan(elevID)));
    if numel(uElevs) == 1
        currentElev = uElevs(1);
        elevStr = sprintf('Elev=%.2f', currentElev);
    elseif numel(uElevs) > 1
        currentElev = mode(uElevs); 
        elevStr = sprintf('Elev=[%.2f..%.2f]', min(uElevs), max(uElevs));
    else
        currentElev = 0;
        elevStr = 'Elev=?';
    end
    
    if isfield(COL,'nHue')
        nHue = COL.nHue;
    else
        nHue = max(hueID(~isnan(hueID)));
    end
    
    hue_all = U.trials.hueID;
    sat_all = U.trials.satID;
    
    saturIDs = unique(satID(~isnan(satID)));
    saturIDs = saturIDs(:).';
    nSat = numel(saturIDs);
    
    if isempty(saturIDs), maxSat = 1; else, maxSat = max(saturIDs); end

    % --- Pre-Calculate Stats ---
    N_hs = zeros(nHue, nSat);
    R_hs = nan(nHue, nSat);
    R_h  = nan(nHue, 1);
    
    for h = 1:nHue
        sel = (hue_all == h);
        R_h(h) = mean(spk_early_norm(sel), 'omitnan');
    end
    
    for si = 1:nSat
        sVal = saturIDs(si);
        selS = (sat_all == sVal);
        for h = 1:nHue
            sel  = selS & (hue_all == h);
            R_hs(h,si) = mean(spk_early_norm(sel), 'omitnan');
            N_hs(h,si) = sum(sel);
        end
    end
    R_h(~isfinite(R_h)) = 0;
    R_hs(~isfinite(R_hs)) = 0;
    
    % --- 4. FIGURE CREATION ---
    figAcc = figure('Name','Spike accounting','Color','w','Visible',figVis);
    figAcc.Units = 'normalized';
    figAcc.Position = [0.1 0.1 0.8 0.8];
    
    nColsTop = 3; 
    nCols    = 3;
    nRowsSat = max(1, ceil(nSat / nCols));
    totalRows = 1 + nRowsSat;
    
    t = tiledlayout(figAcc, totalRows, nCols, 'TileSpacing','normal', 'Padding','compact');
    
    if ~makePlots, set(figAcc, 'Visible', 'off'); end
    
    % --- AXIS LABELING ---
    xtick_vals = 1:nHue;
    xtick_labs = string(xtick_vals);
    xlab_str   = 'Hue Index';
    
    isDKL = false;
    ZERO_HUE = 1; 
    
    if isfield(config,'space') && isfield(config.space,'mode') && strcmpi(config.space.mode,'dkl')
        isDKL = true;
        if isfield(config.space, 'zeroHue'), ZERO_HUE = config.space.zeroHue; end
        
        stepDeg = 360 / nHue;
        % Calculate angle for each hue index relative to Zero Hue
        degrees = mod((xtick_vals - ZERO_HUE) * stepDeg, 360);
        
        % Sparse ticks: Label every 4th hue (assuming 16 total -> 0, 90, 180, 270)
        xtick_labs = strings(1, nHue);
        for i = 1:nHue
            if mod(i-1, 4) == 0 % e.g. indices 1, 5, 9, 13
                d = degrees(i);
                xtick_labs(i) = sprintf('%.0f°', d);
            end
        end
        xlab_str = 'Angle (deg)';
    end
    
    % --- TOP ROW ---
    ax1 = nexttile(t, 1);
    ax2 = nexttile(t, 2);
    ax3 = nexttile(t, 3);
    
    % Panel 1: GROUPED Mean Rates
    axes(ax1); hold(ax1,'on'); %#ok<LAXES>
    b1 = bar(ax1, 1:nHue, R_hs, 'grouped'); 
    for si = 1:nSat
        if si <= numel(b1)
            thisSat = saturIDs(si);
            cDataGroup = get_color_matrix_synthetic(COL, nHue, thisSat, maxSat, currentElev);
            b1(si).CData = cDataGroup;
            b1(si).FaceColor = 'flat';
        end
    end
    title(ax1,'Rate per Hue (Grouped)'); ylabel(ax1,'Norm. Rate');
    grid(ax1,'on'); xlim(ax1,[0.5, nHue+0.5]);
    xticks(ax1, xtick_vals); xticklabels(ax1, xtick_labs); xlabel(ax1, xlab_str);
    
    % Panel 2: GROUPED Trial Counts
    axes(ax2); hold(ax2,'on'); %#ok<LAXES>
    b2 = bar(ax2, 1:nHue, N_hs, 'grouped');
    for si = 1:nSat
        if si <= numel(b2)
            thisSat = saturIDs(si);
            cDataGroup = get_color_matrix_synthetic(COL, nHue, thisSat, maxSat, currentElev);
            b2(si).CData = cDataGroup;
            b2(si).FaceColor = 'flat';
        end
    end
    title(ax2,'Trials per Hue (Grouped)'); ylabel(ax2,'Trials');
    grid(ax2,'on'); xlim(ax2,[0.5, nHue+0.5]);
    xticks(ax2, xtick_vals); xticklabels(ax2, xtick_labs); xlabel(ax2, xlab_str);
    
    % Panel 3: Overall Tuning
    axes(ax3); hold(ax3,'on'); %#ok<LAXES>
    b3 = bar(ax3, 1:nHue, R_h);
    cDataOverall = get_color_matrix_synthetic(COL, nHue, maxSat, maxSat, currentElev);
    b3.CData = cDataOverall;
    b3.FaceColor = 'flat';
    title(ax3,'Overall Tuning'); ylabel(ax3,'Norm. Rate');
    grid(ax3,'on'); xlim(ax3,[0.5, nHue+0.5]);
    xticks(ax3, xtick_vals); xticklabels(ax3, xtick_labs); xlabel(ax3, xlab_str);
    
    % Fit Cosine for Overall + Stats (Amp, Peak, R2)
    hx = linspace(1, nHue, 300);
    x_hue = (0:(nHue-1)).';
    if nHue >= 3
        y_hue = R_h(:);
        Xcos = [ones(nHue,1), cos(2*pi*x_hue/nHue), sin(2*pi*x_hue/nHue)];
        beta = Xcos \ y_hue;
        th_smooth = 2*pi*(hx-1)/nHue; 
        y_recon = beta(1) + beta(2)*cos(th_smooth) + beta(3)*sin(th_smooth);
        plot(ax3, hx, y_recon, 'k-','LineWidth',1.5);
        
        % Calculate Stats
        y_hat = Xcos * beta;
        SS_res = sum((y_hue - y_hat).^2, 'omitnan');
        SS_tot = sum((y_hue - mean(y_hue, 'omitnan')).^2, 'omitnan');
        R2 = 1 - (SS_res / SS_tot);
        if SS_tot < 1e-9, R2 = NaN; end
        
        B = beta(2) + 1i*beta(3);
        
        if isDKL
             % Calculate Peak in DKL Degrees
             pkPh = angle(B); 
             if pkPh < 0, pkPh = pkPh + 2*pi; end
             % Convert Phase (0..2pi) to Index Scale (0..360) then shift by ZeroHue
             pkDegIndex = rad2deg(pkPh); 
             % Shift: If ZeroHue is 16, Index 16 is 0 deg.
             % Angle = (Index - ZeroHue) * Step
             % Phase corresponds to "Index - 1".
             % So: DKL = PhaseDeg - (ZeroHue-1)*Step
             peakDegDisplay = mod(pkDegIndex - (ZERO_HUE-1)*(360/nHue), 360);
        else
             peakDeg = rad2deg(angle(B)); if peakDeg<0, peakDeg=peakDeg+360; end
             peakDegDisplay = peakDeg;
        end
        
        txtStats = sprintf('amp=%.2f, pk=%.0f°, R^2=%.2f', abs(B), peakDegDisplay, R2);
        text(ax3, 0.02, 0.98, txtStats, 'Units','normalized','VerticalAlignment','top','FontSize',9,'EdgeColor','none');
    end
    
    % --- BOTTOM ROWS: Per-Sat Slices ---
    for si = 1:nSat
        tileIdx = nCols + si;
        ax = nexttile(t, tileIdx);
        axes(ax); hold(ax,'on'); %#ok<LAXES>
        
        sVal = saturIDs(si);
        yy = R_hs(:,si);
        bSat = bar(ax, 1:nHue, yy);
        
        cDataSat = get_color_matrix_synthetic(COL, nHue, sVal, maxSat, currentElev);
        bSat.CData = cDataSat;
        bSat.FaceColor = 'flat';
        
        if nHue >= 3
            y_s = yy(:);
            beta_i = Xcos \ y_s;
            y_recon_i = beta_i(1) + beta_i(2)*cos(th_smooth) + beta_i(3)*sin(th_smooth);
            plot(ax, hx, y_recon_i, 'k-','LineWidth',1.0);
            
            % Calc R2 & Stats per sat
            y_hat_s = Xcos * beta_i;
            SS_res_s = sum((y_s - y_hat_s).^2, 'omitnan');
            SS_tot_s = sum((y_s - mean(y_s, 'omitnan')).^2, 'omitnan');
            R2_s = 1 - (SS_res_s / SS_tot_s);
            if SS_tot_s < 1e-9, R2_s = NaN; end
            
            Bi = beta_i(2) + 1i*beta_i(3);
            
            if isDKL
                 pkPh = angle(Bi); if pkPh<0, pkPh=pkPh+2*pi; end
                 pkDegIndex = rad2deg(pkPh);
                 pkDegDisplay = mod(pkDegIndex - (ZERO_HUE-1)*(360/nHue), 360);
            else
                 pkDegI = rad2deg(angle(Bi)); if pkDegI<0, pkDegI=pkDegI+360; end
                 peakDegDisplay = pkDegI;
            end
            
            txtI = sprintf('amp=%.2f, pk=%.0f°, R^2=%.2f', abs(Bi), peakDegDisplay, R2_s);
            text(ax, 0.02, 0.98, txtI, 'Units','normalized','VerticalAlignment','top','FontSize',8,'EdgeColor','none');
        end
        
        title(ax, sprintf('Sat=%.2f', sVal));
        grid(ax,'on'); xlim(ax,[0.5, nHue+0.5]);
        xticks(ax, xtick_vals); xticklabels(ax, xtick_labs); 
        if si > nSat - 3, xlabel(ax, xlab_str); end
    end
    
    % --- TITLE & SAVING ---
    if isfield(U, 'sessionID') && ~isempty(U.sessionID), sessStr = U.sessionID; else, sessStr = dateStr; end
    if isfield(U, 'idx'), idxStr = num2str(U.idx); else, idxStr = '?'; end
    ttl = sprintf('%s | UnitID %d (%s) - UnitNUM %s | %s', sessStr, U.unitID, U.unitType, idxStr, elevStr);
    sgtitle(t, ttl, 'FontWeight','bold', 'Interpreter', 'none');
    
    if ~makePlots, set(figAcc, 'Visible', 'off'); end
    
    baseName = sprintf('%s_%s_%d_spikeaccounting', dateStr, unitType, unitID);
    exportgraphics(figAcc, fullfile(unitDir, [baseName '.png']), 'Resolution', pngDpi);
    if doSaveFig, savefig(figAcc, fullfile(unitDir, [baseName '.fig'])); end
    
    if ~isempty(sessionDir)
        exportgraphics(figAcc, fullfile(sessionDir, [baseName '.png']), 'Resolution', pngDpi);
        if doSaveFig, savefig(figAcc, fullfile(sessionDir, [baseName '.fig'])); end
    end
    if ~isempty(globalDir)
        exportgraphics(figAcc, fullfile(globalDir, [baseName '.png']), 'Resolution', pngDpi);
        if doSaveFig, savefig(figAcc, fullfile(globalDir, [baseName '.fig'])); end
    end
    
    if ~makePlots, close(figAcc); end
end

% -------------------------------------------------------------------------
% Helper: Synthetic Color Scaling (Guarantees consistency)
% -------------------------------------------------------------------------
function cMat = get_color_matrix_synthetic(COL, nHue, sVal, maxSat, eVal)
    cMat = repmat([0.5 0.5 0.5], nHue, 1);
    
    for h = 1:nHue
        % 1. Get Vivid Color (Max Saturation)
        cVivid = get_vivid_color(COL, h, maxSat, eVal);
        
        % 2. Scale brightness by Saturation factor
        % (Assuming Sat=0 is dark/black context, common for bar plots on white)
        if maxSat > 0
            scaleFactor = sVal / maxSat;
        else
            scaleFactor = 0;
        end
        
        cMat(h, :) = cVivid * scaleFactor; 
    end
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