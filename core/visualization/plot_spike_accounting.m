function plot_spike_accounting(U, stats, COL, config, outDir, varargin)
% PLOT_SPIKE_ACCOUNTING - Visual QC of spike counts per condition.
% FIXES: 
% 1. Defines 'doSaveFig' variable (prevents crash).
% 2. Double-locks visibility to 'off' (prevents pop-ups).
% 3. Prioritizes RGB values directly from U.trials.

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
    
    % --- FIX: DEFINE doSaveFig HERE ---
    doSaveFig = true;
    if isfield(config, 'plot') && isfield(config.plot, 'saveFigs')
        doSaveFig = config.plot.saveFigs;
    end
    % ----------------------------------

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
    
    % --- Pre-Calculate Stats ---
    N_hs = zeros(nHue, nSat);
    R_hs = nan(nHue, nSat);
    R_h  = nan(nHue, 1);
    
    % Overall
    for h = 1:nHue
        sel = (hue_all == h);
        R_h(h) = mean(spk_early_norm(sel), 'omitnan');
    end
    
    % Per Sat
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
    
    t = tiledlayout(figAcc, totalRows, nCols, 'TileSpacing','compact', 'Padding','compact');
    
    if ~makePlots, set(figAcc, 'Visible', 'off'); end
    
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
            cDataGroup = get_color_matrix(U.trials, COL, nHue, thisSat, currentElev);
            b1(si).CData = cDataGroup;
            b1(si).FaceColor = 'flat';
        end
    end
    xlabel(ax1,'Hue'); ylabel(ax1,'Norm. Rate'); title(ax1,'Rate per Hue (Grouped)');
    grid(ax1,'on'); xlim(ax1,[0.5, nHue+0.5]);
    
    % Panel 2: GROUPED Trial Counts
    axes(ax2); hold(ax2,'on'); %#ok<LAXES>
    b2 = bar(ax2, 1:nHue, N_hs, 'grouped');
    for si = 1:nSat
        if si <= numel(b2)
            thisSat = saturIDs(si);
            cDataGroup = get_color_matrix(U.trials, COL, nHue, thisSat, currentElev);
            b2(si).CData = cDataGroup;
            b2(si).FaceColor = 'flat';
        end
    end
    xlabel(ax2,'Hue'); ylabel(ax2,'Trials'); title(ax2,'Trials per Hue (Grouped)');
    grid(ax2,'on'); xlim(ax2,[0.5, nHue+0.5]);
    
    % Panel 3: Overall Tuning
    axes(ax3); hold(ax3,'on'); %#ok<LAXES>
    b3 = bar(ax3, 1:nHue, R_h);
    
    % For overall, we use Max Saturation colors (usually mostly vibrant)
    maxSat = max(saturIDs);
    cDataOverall = get_color_matrix(U.trials, COL, nHue, maxSat, currentElev);
    
    b3.CData = cDataOverall;
    b3.FaceColor = 'flat';
    xlabel(ax3,'Hue'); ylabel(ax3,'Norm. Rate'); title(ax3,'Overall Tuning');
    grid(ax3,'on'); xlim(ax3,[0.5, nHue+0.5]);
    
    % Fit Cosine for Overall
    hx = linspace(1, nHue, 300);
    x_hue = (0:(nHue-1)).';
    if nHue >= 3
        y_hue = R_h(:);
        Xcos = [ones(nHue,1), cos(2*pi*x_hue/nHue), sin(2*pi*x_hue/nHue)];
        beta = Xcos \ y_hue;
        th_smooth = 2*pi*(hx-1)/nHue; 
        y_recon = beta(1) + beta(2)*cos(th_smooth) + beta(3)*sin(th_smooth);
        plot(ax3, hx, y_recon, 'k-','LineWidth',1.5);
        
        % Stats
        B = beta(2) + 1i*beta(3);
        peakDeg = rad2deg(angle(B)); if peakDeg < 0, peakDeg=peakDeg+360; end
        txtStats = sprintf('amp=%.2f, peak=%d°', abs(B), round(peakDeg));
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
        
        % Get colors specifically for this Saturation
        cDataSat = get_color_matrix(U.trials, COL, nHue, sVal, currentElev);
        
        bSat.CData = cDataSat;
        bSat.FaceColor = 'flat';
        
        if nHue >= 3
            y_s = yy(:);
            beta_i = Xcos \ y_s;
            y_recon_i = beta_i(1) + beta_i(2)*cos(th_smooth) + beta_i(3)*sin(th_smooth);
            plot(ax, hx, y_recon_i, 'k-','LineWidth',1.0);
            
            Bi = beta_i(2) + 1i*beta_i(3);
            pkDegI = rad2deg(angle(Bi)); if pkDegI<0, pkDegI=pkDegI+360; end
            txtI = sprintf('amp=%.2f, pk=%d°', abs(Bi), round(pkDegI));
            text(ax, 0.02, 0.98, txtI, 'Units','normalized','VerticalAlignment','top','FontSize',8,'EdgeColor','none');
        end
        
        xlabel(ax,'Hue'); 
        title(ax, sprintf('Sat=%.2f', sVal));
        grid(ax,'on'); xlim(ax,[0.5, nHue+0.5]);
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
% Helper: Get entire color matrix (nHue x 3) for a given Saturation/Elev
% Priority: U.trials table > COL metadata
% -------------------------------------------------------------------------
function cMat = get_color_matrix(trials, COL, nHue, sVal, eVal)
    cMat = repmat([0.5 0.5 0.5], nHue, 1);
    
    % Check if trials table has color info (R,G,B or rgb)
    hasRGB = false;
    if any(ismember(trials.Properties.VariableNames, {'R','G','B'}))
        hasRGB = true;
        cMode = 'Split'; % R, G, B columns
    elseif any(ismember(trials.Properties.VariableNames, {'rgb','RGB','color','Color'}))
        hasRGB = true;
        cMode = 'Vector'; % Single column with [r g b] vectors
    end
    
    for h = 1:nHue
        % Default fallback
        c = [0.5 0.5 0.5];
        
        if hasRGB
            % 1. Try fetching from Trials directly
            % Find trials matching this H, S, E
            mask = (trials.hueID == h) & (abs(trials.satID - sVal) < 1e-4);
            
            if any(mask)
                firstIdx = find(mask, 1);
                if strcmp(cMode, 'Split')
                    c = [trials.R(firstIdx), trials.G(firstIdx), trials.B(firstIdx)];
                else
                    % Handle different vector names
                    if ismember('rgb', trials.Properties.VariableNames)
                        c = trials.rgb(firstIdx, :);
                    elseif ismember('RGB', trials.Properties.VariableNames)
                        c = trials.RGB(firstIdx, :);
                    elseif ismember('color', trials.Properties.VariableNames)
                        c = trials.color(firstIdx, :);
                    end
                end
                
                % Normalize if 0-255
                if max(c) > 1.05
                    c = c / 255;
                end
            else
                % If no trials found for this specific combo, fallback to meta
                c = get_color_from_meta_robust(COL, h, sVal, eVal);
            end
        else
            % 2. No RGB in table? Use Robust Metadata Lookup
            c = get_color_from_meta_robust(COL, h, sVal, eVal);
        end
        
        cMat(h, :) = c;
    end
end

function c = get_color_from_meta_robust(COL, h, s, e)
    c = [0.5 0.5 0.5];
    if ~isfield(COL, 'probeIDs') || ~isfield(COL, 'probeCols'), return; end
    tol = 1e-4;
    
    % 1. Exact
    idx = find(abs(COL.probeIDs(:,1)-h)<tol & abs(COL.probeIDs(:,2)-s)<tol & abs(COL.probeIDs(:,3)-e)<tol, 1);
    
    % 2. Fallback (Ignore Elev - treat as Equator)
    if isempty(idx)
        idx = find(abs(COL.probeIDs(:,1)-h)<tol & abs(COL.probeIDs(:,2)-s)<tol & abs(COL.probeIDs(:,3)-0)<tol, 1);
    end
    
    % 3. Desperate (Ignore Elev & Sat match - just find Hue? Risk of wrong saturation color)
    if isempty(idx)
         idx = find(abs(COL.probeIDs(:,1)-h)<tol & abs(COL.probeIDs(:,2)-s)<tol, 1);
    end

    if ~isempty(idx)
        rawC = double(COL.probeCols(idx,:));
        if max(rawC) > 1, rawC = rawC / 255; end
        c = rawC;
    end
end