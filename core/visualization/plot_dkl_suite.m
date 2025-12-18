function plot_dkl_suite(U, stats, COL, config, outDir)
% PLOT_DKL_SUITE - 3D/2D visualization of tuning in DKL space.
% Compatible with sphere_slicer (Robust Lookup + Flat Saving).

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

    % Session Dir (Smart Flat Saving)
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

    % Global Dir
    globalDir = '';
    if isfield(config,'paths') && isfield(config.paths,'globalDiscProbeFigRoot') ...
            && ~isempty(config.paths.globalDiscProbeFigRoot)
        globalDir = fullfile(config.paths.globalDiscProbeFigRoot, 'dkl_suite');
        if ~exist(globalDir,'dir'), mkdir(globalDir); end
    end

    % --- DATA AGGREGATION ---
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
    
    uH = unique(hueID(~isnan(hueID)));
    uS = unique(satID(~isnan(satID)));
    uE = unique(elevID(~isnan(elevID)));
    if isempty(uE), defE = 0; else, defE = mode(uE); end

    X = []; Y = []; Z = []; R = []; C = [];
    r_per_hue = zeros(numel(uH), 1); 
    
    for hi = 1:numel(uH)
        hVal = uH(hi);
        hueRates = [];
        for si = 1:numel(uS)
            sVal = uS(si);
            mask = (hueID == hVal) & (satID == sVal);
            if nnz(mask) == 0, continue; end
            
            meanR = mean(rates(mask), 'omitnan');
            hueRates = [hueRates; meanR]; %#ok<AGROW>
            
            % Robust Lookup
            [dklCoords, rgbColor] = get_coords_from_meta(COL, hVal, sVal, defE);
            if all(isnan(dklCoords)), continue; end
            
            X = [X; dklCoords(1)]; %#ok<AGROW>
            Y = [Y; dklCoords(2)]; %#ok<AGROW>
            Z = [Z; dklCoords(3)]; %#ok<AGROW>
            R = [R; meanR];        %#ok<AGROW>
            C = [C; rgbColor];     %#ok<AGROW>
        end
        if ~isempty(hueRates), r_per_hue(hi) = mean(hueRates); end
    end
    
    if isempty(R), return; end
    C = min(max(C,0),1);
    
    % --- PROJECTION ---
    coords = [X, Y, Z];
    stds   = std(coords, 0, 1);
    [~, minDim] = min(stds); 
    dims = setdiff(1:3, minDim);
    ixX = dims(1); ixY = dims(2);
    XY = coords(:, [ixX, ixY]);

    % --- PLOTTING ---
    f = figure('Color','w','Visible',figVis,'Name','DKL suite');
    f.Units = 'normalized'; f.Position = [0.1 0.1 0.75 0.7];
    tlo = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');

    % 3D Bars
    ax3d = nexttile(tlo,1);
    hold(ax3d,'on'); grid(ax3d,'on'); axis(ax3d,'equal');
    rPos = max(R,0); maxR = max(rPos); if maxR <= 0, rPos(:) = 0; end
    baseRadius = 0.05 * (max(XY(:)) - min(XY(:)));
    if baseRadius == 0, baseRadius = 0.05; end
    
    for k = 1:numel(rPos)
        if rPos(k) <= 0, continue; end
        x0 = XY(k,1); y0 = XY(k,2); z0 = 0; h = rPos(k); w = baseRadius;
        xv = [x0-w x0+w x0+w x0-w]; yv = [y0-w y0-w y0+w y0+w];
        zv = [z0   z0   z0   z0];   zt = zv + h;
        fc = C(k,:);
        patch(ax3d, xv, yv, zt, fc, 'FaceAlpha',0.9, 'EdgeColor','none');
    end
    xlabel(ax3d, sprintf('DKL Axis %d', ixX)); ylabel(ax3d, sprintf('DKL Axis %d', ixY));
    zlabel(ax3d, 'Rate (Hz)'); title(ax3d, 'Response in DKL Plane'); view(ax3d, [40 25]);

    % 2D Summary
    ax2d = nexttile(tlo,2); hold(ax2d,'on'); grid(ax2d,'on');
    plot(ax2d, 1:numel(uH), r_per_hue, 'o-', 'LineWidth', 1.5, 'Color', 'k');
    xlim(ax2d, [0.5 numel(uH)+0.5]);
    xlabel(ax2d, 'Hue Index'); ylabel(ax2d, 'Mean Rate (Hz)'); title(ax2d, 'Marginal Hue Tuning');

    % Title Handling (Replaced dash with standard hyphen)
    if isfield(U,'exptName')
        ttl = sprintf('%s | Unit %d (%s) - DKL suite', string(U.exptName), U.unitID, U.unitType);
    else
        ttl = sprintf('%s | Unit %d (%s) - DKL suite', string(U.dateStr), U.unitID, U.unitType);
    end
    sgtitle(f, ttl, 'Interpreter','none');

    % --- SAVING ---
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
    if ~isempty(idx)
        dkl = double(COL.dklRows(idx, :));
        if isfield(COL, 'probeCols'), rawC = double(COL.probeCols(idx, :)); if max(rawC)>1, rawC=rawC/255; end; rgb=rawC; end
    end
end