function plot_dkl_suite(U, stats, COL, config, outDir)
% PLOT_DKL_SUITE - Polar visualization of tuning in DKL space.
% FINAL FIXES:
% 1. Zero Alignment: Uses config.space.zeroHue to set which Hue is at 0°.
% 2. Spacing: Pushes Labels OUTSIDE the Hue Ring.
% 3. Cleanliness: Removes redundant subplot title.

    if U.nTrials == 0
        return;
    end
    
    % --- CONFIG ---
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
            sessionDir = fullfile(outDir, 'dkl_suite');
        else
            sessionDir = outDir; 
        end
        if ~exist(sessionDir,'dir'), mkdir(sessionDir); end
    end

    globalDir = '';
    if isfield(config,'paths') && isfield(config.paths,'globalDiscProbeFigRoot') ...
            && ~isempty(config.paths.globalDiscProbeFigRoot)
        globalDir = fullfile(config.paths.globalDiscProbeFigRoot, 'dkl_suite');
        if ~exist(globalDir,'dir'), mkdir(globalDir); end
    end

    % --- DATA PREP ---
    win = U.winEarly;
    dur = diff(win);
    
    nTr = numel(U.spk);
    rates = nan(nTr,1);
    for tt = 1:nTr
        spks = U.spk{tt};
        if isempty(spks)
            rates(tt) = 0;
        else
            spks = spks(spks >= win(1) & spks < win(2));
            rates(tt) = numel(spks) / dur;
        end
    end
    
    hueID = U.trials.hueID; satID = U.trials.satID; elevID = U.trials.elevID;
    
    if isfield(COL, 'nHue')
        nHue = COL.nHue;
    else
        nHue = max(hueID(~isnan(hueID)));
    end
    
    uS = unique(satID(~isnan(satID)));
    uE = unique(elevID(~isnan(elevID)));
    if isempty(uE), defE = 0; else, defE = mode(uE); end
    if isempty(uS), maxS = 1; else, maxS = max(uS); end

    % --- ROTATION FIX: Handle Zero Hue ---
    ZERO_HUE = 1; % Default
    if isfield(config, 'space') && isfield(config.space, 'zeroHue')
        ZERO_HUE = config.space.zeroHue;
    end

    hue_indices = (1:nHue)';
    
    % Calculate Angle relative to Zero Hue
    % mod(h - zero, n) ensures correct rotation
    % e.g. if Zero=16, Hue 16 becomes 0 deg, Hue 1 becomes 22.5 deg
    hue_angles = 2*pi * mod(hue_indices - ZERO_HUE, nHue) / nHue;
    
    r_per_hue = zeros(nHue, 1);
    for hVal = 1:nHue
        mask = (hueID == hVal);
        if nnz(mask) > 0
            r_per_hue(hVal) = mean(rates(mask), 'omitnan');
        else
            r_per_hue(hVal) = 0;
        end
    end
    
    % --- PLOTTING ---
    f = figure('Color','w','Visible',figVis,'Name','DKL suite');
    f.Units = 'normalized'; f.Position = [0.2 0.2 0.5 0.75]; 
    
    axP = polaraxes('Parent', f);
    hold(axP, 'on');
    
    % 1. ORIENTATION
    axP.ThetaZeroLocation = 'right';
    axP.ThetaDir = 'counterclockwise';
    axP.FontSize = 10; 
    
    % Sort by angle for proper line drawing
    [sortedAngles, sortIdx] = sort(hue_angles);
    sortedRates = r_per_hue(sortIdx);
    sortedHues  = hue_indices(sortIdx);
    
    % Close loop
    if ~isempty(sortedAngles)
        sortedAngles(end+1) = sortedAngles(1) + 2*pi;
        sortedRates(end+1)  = sortedRates(1);
    end
    
    % --- SCALING LAYERS ---
    maxRate = max(sortedRates, [], 'omitnan');
    if maxRate <= 0, maxRate = 1; end
    
    ringCenter  = maxRate * 1.15;  
    viewLimit   = maxRate * 1.25;  
    labelRadius = maxRate * 1.55; 
    
    rlim(axP, [0 viewLimit]);
    
    % 2. Draw Tuning Curve (Red)
    polarplot(axP, sortedAngles, sortedRates, 'o-', ...
        'LineWidth', 2.5, 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
    
    % 3. Draw Colored Ring (Flush with edge)
    hueWidth = 2*pi / nHue;
    
    for k = 1:nHue
        % Use sorted indices to match the angles
        hVal = sortedHues(k);
        thisAngle = sortedAngles(k);
        
        [~, ringC] = get_coords_from_meta(COL, hVal, maxS, defE);
        
        thStart = thisAngle - hueWidth/2 * 0.98;
        thEnd   = thisAngle + hueWidth/2 * 0.98;
        thArc   = linspace(thStart, thEnd, 20); 
        rArc    = repmat(ringCenter, size(thArc));
        
        polarplot(axP, thArc, rArc, '-', 'Color', ringC, 'LineWidth', 20); 
    end
    
    % 4. Preferred Direction
    if ~isempty(r_per_hue)
        % Calculate vector sum using original angles
        V = sum(r_per_hue .* exp(1i * hue_angles));
        prefAng = angle(V);
        polarplot(axP, [0, prefAng], [0, maxRate], '--', ...
            'Color', [0.9 0.7 0], 'LineWidth', 2);
    end
    
    % 5. DKL Labels (Outside everything)
    fontSz = 11;
    text(axP, 0,    labelRadius, '+[L-M]', 'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize', fontSz, 'Clipping', 'off');
    text(axP, pi/2, labelRadius, '+S',     'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize', fontSz, 'Clipping', 'off');
    text(axP, pi,   labelRadius, '-[L-M]', 'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize', fontSz, 'Clipping', 'off');
    text(axP, 3*pi/2, labelRadius, '-S',   'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize', fontSz, 'Clipping', 'off');

    % --- TITLE & SAVING ---
    if isfield(U,'exptName')
        ttl = sprintf('%s | Unit %d (%s) - DKL Polar', string(U.exptName), U.unitID, U.unitType);
    else
        ttl = sprintf('%s | Unit %d (%s) - DKL Polar', string(U.dateStr), U.unitID, U.unitType);
    end
    
    % Safe title positioning (axes title with newlines)
    title(axP, {ttl, ''}, 'Interpreter','none', 'FontSize', 12, 'FontWeight', 'bold');

    fileTag = sprintf('%s_%s_%d', string(U.dateStr), U.unitType, U.unitID);
    exportgraphics(f, fullfile(unitDir, sprintf('DKL_suite_%s.png', fileTag)), 'Resolution', pngDpi);
    if doSaveFig, savefig(f, fullfile(unitDir, sprintf('DKL_suite_%s.fig', fileTag))); end
    
    if ~isempty(sessionDir)
        exportgraphics(f, fullfile(sessionDir, sprintf('DKL_suite_%s.png', fileTag)), 'Resolution', pngDpi);
        if doSaveFig, savefig(f, fullfile(sessionDir, sprintf('DKL_suite_%s.fig', fileTag))); end
    end
    
    if ~isempty(globalDir)
        exportgraphics(f, fullfile(globalDir, sprintf('DKL_suite_%s.png', fileTag)), 'Resolution', pngDpi);
        if doSaveFig, savefig(f, fullfile(globalDir, sprintf('DKL_suite_%s.fig', fileTag))); end
    end
    
    if ~makePlots, close(f); end
end

% --- HELPER ---
function [dkl, rgb] = get_coords_from_meta(COL, h, s, e)
    dkl = [NaN NaN NaN]; rgb = [0.5 0.5 0.5];
    if ~isfield(COL, 'probeIDs') || ~isfield(COL, 'dklRows'), return; end
    tol = 1e-4;
    idx = find(abs(COL.probeIDs(:,1)-h)<tol & abs(COL.probeIDs(:,2)-s)<tol & abs(COL.probeIDs(:,3)-e)<tol, 1);
    
    % Fallback 1: Ignore Elev
    if isempty(idx)
         idx = find(abs(COL.probeIDs(:,1)-h)<tol & abs(COL.probeIDs(:,2)-s)<tol, 1);
    end
    % Fallback 2: Any Hue match (useful for Vivid lookup)
    if isempty(idx)
        idx = find(abs(COL.probeIDs(:,1)-h)<tol, 1);
    end
    
    if ~isempty(idx)
        dkl = double(COL.dklRows(idx, :));
        if isfield(COL, 'probeCols'), rawC = double(COL.probeCols(idx, :)); if max(rawC)>1, rawC=rawC/255; end; rgb=rawC; end
    end
end