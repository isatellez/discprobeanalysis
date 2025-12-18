function H = plot_pooled_wachtler_fits_histograms(config)
% Pool all *_wachtler_fits.csv files and plot population histograms

    if nargin < 1 || isempty(config)
        config = load_config();
    end

    H = struct('fig', [], 'axDir', [], 'axWidth', [], ...
               'pngPath', [], 'pdfPath', []);

    if isfield(config,'paths') && isfield(config.paths,'output') && ~isempty(config.paths.output)
        outRoot = config.paths.output;
    else
        rootDir = fileparts(mfilename('fullpath'));
        outRoot = fullfile(rootDir, 'output');
    end

    wachtRoot = fullfile(outRoot, 'wachtler');
    if ~isfolder(wachtRoot)
        error('plot_pooled_wachtler_fits_histograms:NoWachtlerFolder', ...
            'No wachtler folder found at %s', wachtRoot);
    end

    D = dir(fullfile(wachtRoot, '**', '*_wachtler_fits.csv'));
    if isempty(D)
        error('plot_pooled_wachtler_fits_histograms:NoCSV', ...
            'No *_wachtler_fits.csv files found under %s', wachtRoot);
    end

    allPhi  = [];
    allFwhm = [];
    nFilesUsed = 0;

    for k = 1:numel(D)
        [dateFolder, ~] = fileparts(D(k).folder);
        [~, dateName]   = fileparts(dateFolder);

        tok = regexp(dateName, '^\d{6}$', 'match', 'once');
        if isempty(tok)
            continue;
        end

        csvPath = fullfile(D(k).folder, D(k).name);
        try
            T = readtable(csvPath);
        catch ME
            warning('Could not read %s: %s', csvPath, ME.message);
            continue;
        end

        vnames = T.Properties.VariableNames;
        needed = ismember({'prefDeg','FWHMDeg','R2'}, vnames);
        if ~all(needed)
            warning('Skipping %s (missing required columns).', csvPath);
            continue;
        end

        if ismember('isTunedR2', vnames)
            tuned = T.isTunedR2 ~= 0;
        else
            if ~isfield(config,'tuning') || ~isfield(config.tuning,'R2_thresh') ...
                    || isempty(config.tuning.R2_thresh)
                config.tuning.R2_thresh = 0.5;
            end
            tuned = T.R2 >= config.tuning.R2_thresh;
        end

        tuned = tuned & isfinite(T.prefDeg) & isfinite(T.FWHMDeg);

        allPhi  = [allPhi;  T.prefDeg(tuned)];
        allFwhm = [allFwhm; T.FWHMDeg(tuned)];
        nFilesUsed = nFilesUsed + 1;
    end

    if isempty(allPhi)
        warning('No tuned units found in any wachtler_fits files.');
        return;
    end

    fprintf('Pooled Wachtler fits from %d CSV files (%d tuned units).\n', ...
        nFilesUsed, numel(allPhi));

    % figure roughly like the paper
    H.fig = figure('Color','w', ...
                   'Units','centimeters', ...
                   'Position',[2 2 9 11]);  % [left bottom width height]

    %% ---- Panel A: tuning peak direction ----
    % ---- Panel A: tuning peak direction ----
    H.axDir = subplot(2,1,1);

    allPhiWrapped = mod(allPhi, 360);

    % 16 equal bins, centers at 0, 22.5, ..., 337.5
    nBinsDir     = 16;
    binWidthDir  = 360 / nBinsDir;      % 22.5
    halfWidthDir = binWidthDir / 2;     % 11.25

    % edges: -11.25, 11.25, ..., 348.75  -> 16 bins total
    edgesDir   = -halfWidthDir : binWidthDir : (360 - halfWidthDir);
    centersDir = edgesDir(1:end-1) + halfWidthDir;   % 0:22.5:337.5

    countsDir = histcounts(allPhiWrapped, edgesDir);

    bar(centersDir, countsDir, 0.9, ...   % 0.9 so the y-axis isn't hiding the bar
        'FaceColor', [0.8 0.8 0.8], ...
        'EdgeColor', 'k');

    ymax = max(countsDir) + 1;
    ylim([0 ymax]);
    yticks(0:5:ymax);

    xlim([0 360]);
    xticks(0:90:360);
    xticklabels({'0','90','180','270','360'});

    ylabel('number of cells');
    xlabel('tuning peak direction [deg]');

    set(H.axDir, 'Box','off', 'TickDir','out', 'LineWidth',1);
    set(H.axDir, 'Clipping','off');

    % put +L-M, +S, -L+M, -S in *normalized* coords, well below degree ticks
    hold on;
    yNorm = -0.30;   % more negative = further below x-axis / tick labels

    text(0/360,   yNorm, '+L-M', 'Units','normalized', ...
         'HorizontalAlignment','center');
    text(90/360,  yNorm, '+S',   'Units','normalized', ...
         'HorizontalAlignment','center');
    text(180/360, yNorm, '-L+M', 'Units','normalized', ...
         'HorizontalAlignment','center');
    text(270/360, yNorm, '-S',   'Units','normalized', ...
         'HorizontalAlignment','center');
    hold off;


    %% ---- Panel B: tuning width ----
    H.axWidth = subplot(2,1,2);

    edgesW   = 0:15:180;                  % 0–15, 15–30, ..., 165–180
    centersW = edgesW(1:end-1) + diff(edgesW)/2;

    % keep only finite widths in (0, 180]; discard weird / out-of-range ones
    validFwhm = allFwhm;
    validFwhm = validFwhm(isfinite(validFwhm) & validFwhm > 0 & validFwhm <= 180);

    countsW = histcounts(validFwhm, edgesW);

    bar(centersW, countsW, 1.0, ...
        'FaceColor', [0.8 0.8 0.8], 'EdgeColor','k');

    xlim([0 180]);
    xticks(0:30:180);
    xlabel('tuning width [deg]');
    ylabel('number of cells');

    set(H.axWidth, 'Box','off', 'TickDir','out', 'LineWidth',1);
    hold on;
    xline(120, 'k:');
    hold off;


    % save figure (PNG + PDF) in wachtler root
    figBase = fullfile(wachtRoot, 'wachtler_pooled_histograms');
    H.pngPath = [figBase '.png'];
    H.pdfPath = [figBase '.pdf'];

    set(H.fig, 'PaperUnits','centimeters', 'PaperPosition',[0 0 9 11]);
    print(H.fig, H.pngPath, '-dpng', '-r300');
    print(H.fig, H.pdfPath, '-dpdf', '-r300');
end
