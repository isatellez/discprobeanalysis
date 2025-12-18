function plot_planes_rayleigh_spike(sessionID, config)
% For each isoluminant plane (elevID), plot:
%   - Rayleigh R vs saturation (with significance stars)
%   - best-model R^2 vs saturation (spike accounting)

    if nargin < 2 || isempty(config)
        config = load_config();
    end

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

    % Rayleigh alpha (per-unit)
    alphaRay = 0.05;
    if isfield(config,'tuning') && isfield(config.tuning,'alphaRayleigh') ...
            && ~isempty(config.tuning.alphaRayleigh)
        alphaRay = config.tuning.alphaRayleigh;
    end

    P = get_session_paths(config, sessionID);
    tablesDir = P.tables;
    figsDir   = P.figs;

    sessTag = char(sessionID);
    csvFile = fullfile(tablesDir, sprintf('%s_sphereSlices.csv', sessTag));

    if ~isfile(csvFile)
        warning('plot_planes_rayleigh_spike: CSV not found: %s', csvFile);
        return;
    end

    T = readtable(csvFile);
    if isempty(T)
        warning('plot_planes_rayleigh_spike: empty table in %s', csvFile);
        return;
    end

    if ~ismember('elevID', T.Properties.VariableNames) || ...
       ~ismember('satID',  T.Properties.VariableNames)
        error('plot_planes_rayleigh_spike: CSV missing satID/elevID.');
    end

    % bestR2 column may or may not exist; if not, build it from r2_1/r2_12
    if ismember('bestR2', T.Properties.VariableNames)
        bestR2_all = T.bestR2;
    else
        r2_1_all  = NaN(height(T),1);
        r2_12_all = NaN(height(T),1);
        if ismember('r2_1', T.Properties.VariableNames)
            r2_1_all = T.r2_1;
        end
        if ismember('r2_12', T.Properties.VariableNames)
            r2_12_all = T.r2_12;
        end
        bestR2_all = r2_1_all;
        use12 = isfinite(r2_12_all) & (r2_12_all > r2_1_all | ~isfinite(r2_1_all));
        bestR2_all(use12) = r2_12_all(use12);
        T.bestR2 = bestR2_all;
    end

    elevVals = unique(T.elevID);
    satVals  = unique(T.satID);
    satVals  = satVals(isfinite(satVals));

    outDir = fullfile(figsDir, 'discprobe', 'planes');
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end

    globalDir = '';
    if isfield(config,'paths') && isfield(config.paths,'globalDiscProbeFigRoot') ...
            && ~isempty(config.paths.globalDiscProbeFigRoot)
        globalDir = fullfile(config.paths.globalDiscProbeFigRoot, 'planes');
        if ~exist(globalDir,'dir')
            mkdir(globalDir);
        end
    end

    for ei = 1:numel(elevVals)
        e = elevVals(ei);
        usePlane = T.elevID == e;

        if ~any(usePlane)
            continue;
        end

        satList = satVals;
        mRay        = nan(size(satList));
        seRay       = nan(size(satList));
        mBestR2     = nan(size(satList));
        seBestR2    = nan(size(satList));
        fracTuned   = nan(size(satList));
        nUnits      = zeros(size(satList));
        nUnitsTuned = zeros(size(satList));

        for si = 1:numel(satList)
            s = satList(si);
            use = usePlane & abs(T.satID - s) < 1e-6;

            if ~any(use)
                continue;
            end

            rayVals = T.rayleigh_R(use);
            r2Vals  = T.bestR2(use);

            % Rayleigh p-values for significance
            rayP = NaN(sum(use),1);
            if ismember('rayleigh_p', T.Properties.VariableNames)
                rayP = T.rayleigh_p(use);
            end

            rayVals = rayVals(isfinite(rayVals));
            r2Vals  = r2Vals(isfinite(r2Vals));

            if ~isempty(rayVals)
                mRay(si)  = mean(rayVals);
                seRay(si) = std(rayVals) / sqrt(numel(rayVals));
            end
            if ~isempty(r2Vals)
                mBestR2(si)   = mean(r2Vals);
                seBestR2(si)  = std(r2Vals) / sqrt(numel(r2Vals));
            end

            if ~isempty(rayP)
                tunedMask = isfinite(rayP) & rayP < alphaRay;
                nUnits(si)      = numel(rayP);
                nUnitsTuned(si) = sum(tunedMask);
                if nUnits(si) > 0
                    fracTuned(si) = nUnitsTuned(si) / nUnits(si);
                end
            end
        end

        fig = figure('Name', sprintf('%s elev=%g', sessTag, e), ...
                     'Color','w','Visible',figVis);
        fig.Units = 'normalized';
        fig.Position = [0.2 0.2 0.5 0.6];

        % ---- Rayleigh vs saturation ----
        ax1 = subplot(2,1,1);
        hold(ax1,'on');
        errorbar(ax1, satList, mRay, seRay, 'o-', 'LineWidth', 1.5);
        xlabel(ax1, 'saturation');
        ylabel(ax1, 'Rayleigh R');
        title(ax1, sprintf('%s: Rayleigh vs sat (elev=%g)', sessTag, e), ...
              'Interpreter','none');
        grid(ax1,'on');

        % place stars where there are any tuned units at that sat/plane
        yMax = max(mRay + seRay, [], 'omitnan');
        if isfinite(yMax) && yMax > 0
            yStar = yMax * 1.05;
            for si = 1:numel(satList)
                if nUnitsTuned(si) > 0
                    plot(ax1, satList(si), yStar, 'k*', 'MarkerSize', 6);
                end
            end
            ylim(ax1, [0, yStar * 1.1]);
        end

        % small legend note for alpha
        txtAlpha = sprintf('\\alpha_{Rayleigh} = %.3f; star = any tuned units', alphaRay);
        text(ax1, 0.01, 0.97, txtAlpha, ...
            'Units','normalized', ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','top');

        % ---- best-model R^2 vs saturation ----
        ax2 = subplot(2,1,2);
        hold(ax2,'on');
        errorbar(ax2, satList, mBestR2, seBestR2, 'o-', 'LineWidth', 1.5);
        xlabel(ax2, 'saturation');
        ylabel(ax2, 'best-model R^2');
        title(ax2, sprintf('%s: spike accounting (best R^2) vs sat (elev=%g)', sessTag, e), ...
              'Interpreter','none');
        grid(ax2,'on');

        % optional text with how many units per sat
        for si = 1:numel(satList)
            if nUnits(si) > 0
                txt = sprintf('%d/%d tuned', nUnitsTuned(si), nUnits(si));
                text(ax2, satList(si), mBestR2(si), txt, ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','bottom', ...
                    'FontSize',8);
            end
        end

        sgtitle(fig, sprintf('%s – plane elev=%g', sessTag, e), 'Interpreter','none');

        baseName = sprintf('%s_plane_elev%+0.2f', sessTag, e);
        png_plane = fullfile(outDir, [baseName '.png']);
        fig_plane = fullfile(outDir, [baseName '.fig']);
        exportgraphics(fig, png_plane, 'Resolution', pngDpi);
        savefig(fig, fig_plane);

        if ~isempty(globalDir)
            png_glob = fullfile(globalDir, [baseName '.png']);
            fig_glob = fullfile(globalDir, [baseName '.fig']);
            exportgraphics(fig, png_glob, 'Resolution', pngDpi);
            savefig(fig, fig_glob);
        end

        if ~makePlots
            close(fig);
        end
    end
end
