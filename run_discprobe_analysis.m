function run_discprobe_analysis(DATES)

clc;

rootDir = fileparts(mfilename('fullpath'));
addpath(genpath(rootDir));

config = load_config();
COL    = load_color_metadata(config);

if nargin < 1 || isempty(DATES)
    DATES = {'250513'};
end
DATES = normalize_dates(DATES);

% main output root
if isfield(config, 'paths') && isfield(config.paths, 'output') && ~isempty(config.paths.output)
    outRoot = config.paths.output;
else
    outRoot = fullfile(rootDir, 'output');
end

% RNG setup (unchanged)
if isfield(config, 'pt') && isfield(config.pt, 'seed') && numel(config.pt.seed) == 1
    rng(config.pt.seed);
end

% helper for config.do flags (unchanged)
doFlag = @(name,default) ...
    ( isfield(config,'do') && isfield(config.do,name) && ...
      all(logical(config.do.(name)(:))) ) ...
    || ( ~isfield(config,'do') && default );

for d = 1:numel(DATES)
    dateStr = DATES{d};
    fprintf('\n========== %s ===========\n', dateStr);

    % per-session roots
    sessionRoot = fullfile(outRoot, dateStr);
    figRoot     = fullfile(sessionRoot, 'figs');
    unitsRoot   = fullfile(sessionRoot, 'units');
    tablesRoot  = fullfile(sessionRoot, 'tables');

    if ~isfolder(sessionRoot), mkdir(sessionRoot); end
    if ~isfolder(figRoot),     mkdir(figRoot);     end
    if ~isfolder(unitsRoot),   mkdir(unitsRoot);   end
    if ~isfolder(tablesRoot),  mkdir(tablesRoot);  end

    % if any downstream code still reads config.paths.*, set them here
    config.paths.figRoot   = figRoot;
    config.paths.unitsRoot = unitsRoot;
    config.paths.tablesRoot = tablesRoot;

    % from here down, use figRoot / unitsRoot / tablesRoot
    % instead of the old dateFigRoot/dateDataRoot stuff

    % load session
    S = load_session(dateStr, config);

    % drop units with 0 spikes
    S = filter_zero_spike_units(S, config);

    % trial table (if you want the CSV under tablesRoot, pass it in or
    % change make_trial_index_table to use config.paths.tablesRoot)
    T = make_trial_index_table(S, COL, config);

    Tall = table();

    for unitIdx = 1:S.nUnits
        fprintf('  unit %d/%d (ID=%d)\n', unitIdx, S.nUnits, S.unitIDs(unitIdx));

        U = prepare_unit(S, unitIdx, T, config);

        if doFlag('dropOutliers', false)
            U = detect_outliers(U, config);
        end

        stats = struct();
        stats.hueMeans = compute_hue_means(U);

        if doFlag('permutationTest', false)
            stats.permutation = permutation_test(U, config);
        end
        if doFlag('rayleigh', false)
            stats.rayleigh = rayleigh_test(U, config, COL);
        end
        if doFlag('fourier', false)
            stats.fourier = fourier_analysis(U, config);
        end
        if doFlag('cosineFit', false)
            stats.cosine = cosine_fit(U, config);
        end
        if doFlag('anova', false)
            stats.anova = anova_analysis(U, config);
        end
        if doFlag('peakModel', false)
            stats.peak = peak_model(U, config);
        end

        % plots – now using figRoot
        if doFlag('rasterCanonical', true)
            plot_rasters(U, stats, COL, config, figRoot, "canonical");
        end
        if doFlag('rasterByMean', false)
            plot_rasters(U, stats, COL, config, figRoot, "mean");
        end
        if doFlag('tuningFigure', true)
            plot_tuning_curves(U, stats, COL, config, figRoot);
        end
        if doFlag('rayleigh', false)
            plot_rayleigh(U, stats, COL, config, figRoot);
        end
        if doFlag('fourierFig', false)
            plot_fourier(U, stats, COL, config, figRoot);
        end
        if doFlag('boxPlots', false)
            plot_boxplots(U, stats, COL, config, figRoot);
        end
        if doFlag('spikeAccounting', false)
            plot_spike_accounting(U, stats, COL, config, figRoot);
        end
        if doFlag('dklSuite', true)
            plot_dkl_suite(U, stats, COL, config, figRoot);
        end

        % save per-unit stuff under unitsRoot
        save_unit_outputs(unitsRoot, U, stats);

        % light summary row (same as before)
        if isfield(stats, 'hueMeans') && isfield(stats.hueMeans, 'rate_mean')
            row = table;
            row.dateStr        = string(U.dateStr);
            row.unitID         = U.unitID;
            row.unitIdx        = unitIdx;
            row.unitType       = string(U.unitType);
            row.meanRateOverall = mean(stats.hueMeans.rate_mean, 'omitnan');
            Tall = [Tall; row]; %#ok<AGROW>
        end
    end

    % per-date summary CSV in tablesRoot
    if ~isempty(Tall)
        sumFile = fullfile(tablesRoot, sprintf('%s_unit_summary.csv', dateStr));
        try
            writetable(Tall, sumFile);
        catch ME
            warning('Could not write per-date summary CSV for %s: %s', dateStr, ME.message);
        end
    end
end

fprintf('\nDone.\n');

end
