function run_discprobe_analysis(DATES)
% Main entry point for DiscProbeAnalysis.
% Pulls logic from Analyze_DiscProbePSTHs_v1figproc.m.

clc;

rootDir = fileparts(mfilename('fullpath'));
addpath(genpath(rootDir));

config = load_config();
COL    = load_color_metadata(config);

% dates
if nargin < 1 || isempty(DATES)
    DATES = {'250513'};
end
DATES = normalize_dates(DATES);

% output roots
if isfield(config, 'paths') && isfield(config.paths, 'output') ...
        && ~isempty(config.paths.output)
    outRoot = config.paths.output;
else
    outRoot = fullfile(rootDir, 'output');
end

figRoot  = fullfile(outRoot, 'figs');
dataRoot = fullfile(outRoot, 'summary');

if ~isfolder(figRoot),  mkdir(figRoot);  end
if ~isfolder(dataRoot), mkdir(dataRoot); end

config.paths.figRoot  = figRoot;
config.paths.dataRoot = dataRoot;

% RNG for permutation tests etc
if isfield(config, 'pt') && isfield(config.pt, 'seed') ...
        && numel(config.pt.seed) == 1
    rng(config.pt.seed);
end

% helper to read config.do flags as scalar booleans
doFlag = @(name,default) ...
    ( isfield(config,'do') && isfield(config.do,name) && ...
      all(logical(config.do.(name)(:))) ) ...
    || ( ~isfield(config,'do') && default );

for d = 1:numel(DATES)
    dateStr = DATES{d};
    fprintf('\n========== %s ============\n', dateStr);

    dateFigRoot  = fullfile(figRoot,  dateStr);
    dateDataRoot = fullfile(dataRoot, dateStr);
    if ~isfolder(dateFigRoot),  mkdir(dateFigRoot);  end
    if ~isfolder(dateDataRoot), mkdir(dateDataRoot); end

    %try
        % load session
        S = load_session(dateStr, config);

        % drop units with 0 spikes across full window
        S = filter_zero_spike_units(S, config);

        % trial-level info (RGB → IDs + DKL all handled here)
        T = make_trial_index_table(S, COL, config);

        Tall = table();

        for unitIdx = 1:S.nUnits
            fprintf('  unit %d/%d (ID=%d)\n', ...
                unitIdx, S.nUnits, S.unitIDs(unitIdx));

            % per-unit struct
            U = prepare_unit(S, unitIdx, T, config);

            % mark/drop outliers if requested
            if doFlag('dropOutliers', false)
                U = detect_outliers(U, config);
            end

            % analysis
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

            % visualization
            if doFlag('rasterCanonical', true)
                plot_rasters(U, stats, COL, config, dateFigRoot, "canonical");
            end

            if doFlag('rasterByMean', false)
                plot_rasters(U, stats, COL, config, dateFigRoot, "mean");
            end

            if doFlag('tuningFigure', true)
                plot_tuning_curves(U, stats, COL, config, dateFigRoot);
            end

            if doFlag('rayleigh', false)
                plot_rayleigh(U, stats, COL, config, dateFigRoot);
            end

            if doFlag('boxPlots', false)
                plot_boxplots(U, stats, COL, config, dateFigRoot);
            end

            if doFlag('spikeAccounting', false)
                plot_spike_accounting(U, stats, COL, config, dateFigRoot);
            end

            if doFlag('dklSuite', false)
                plot_dkl_suite(U, stats, COL, config, dateFigRoot);
            end

            % save per-unit outputs
            save_unit_outputs(dateDataRoot, U, stats);

            % light per-date summary row
            if isfield(stats, 'hueMeans') && isfield(stats.hueMeans, 'rate_mean')
                row = table;
                row.dateStr = string(U.dateStr);
                row.unitID  = U.unitID;
                row.unitIdx = unitIdx;
                row.unitType = string(U.unitType);
                row.meanRateOverall = mean(stats.hueMeans.rate_mean, 'omitnan');
                Tall = [Tall; row]; %#ok<AGROW>
            end
        end

        % save per-date summary CSV
        if ~isempty(Tall)
            sumFile = fullfile(dateDataRoot, sprintf('%s_unit_summary.csv', dateStr));
            try
                writetable(Tall, sumFile);
            catch ME
                warning('Could not write per-date summary CSV for %s: %s', dateStr, ME.message);
            end
        end

    % catch ME
    %     warning('Failed on %s: %s', dateStr, ME.message);
    %     fprintf('Skipping %s.\n', dateStr);
    %     continue;
    % end
end


fprintf('\nDone.\n');

end
