function run_discprobe_analysis(DATES)
clc;
set(0,'DefaultFigureVisible','off');
rootDir = fileparts(mfilename('fullpath'));
addpath(genpath(rootDir));

% --- CONFIG ---
config = load_discprobe_config(); 
COL    = load_color_metadata(config);

% default if nothing passed in
if nargin < 1 || isempty(DATES)
    DATES = {'250513'};
end

% handle 'all' and string vs cell input
if ischar(DATES) || isstring(DATES)
    ds = strtrim(char(DATES));
    if strcmpi(ds, 'all')
        if ~isfield(config, 'paths') || ~isfield(config.paths, 'base') || isempty(config.paths.base)
            error('run_discprobe_analysis:NoBasePath', ...
                  'config.paths.base is missing or empty; cannot resolve ''all'' dates.');
        end
        baseRoot = config.paths.base;
        if ~isfolder(baseRoot)
            error('run_discprobe_analysis:BadBasePath', ...
                  'Base path does not exist: %s', baseRoot);
        end
        dInfo = dir(baseRoot);
        names = {dInfo.name};
        
        % Allow Monkey_Date format
        isDateDir = [dInfo.isdir] & ...
                    ~ismember(names, {'.','..'}) & ...
                    ~cellfun(@isempty, regexp(names, '(\d{6})$', 'once'));
                    
        dateNames = names(isDateDir);
        if isempty(dateNames)
            error('run_discprobe_analysis:NoDateFolders', ...
                  'No valid date folders found under %s', baseRoot);
        end
        DATES = sort(dateNames);
    else
        DATES = {ds};
    end
end
DATES = normalize_dates(DATES);

% RNG setup
if isfield(config, 'pt') && isfield(config.pt, 'seed') && numel(config.pt.seed) == 1
    rng(config.pt.seed);
end

% global "overall" discprobe roots
if isfield(config,'paths') && isfield(config.paths,'output') && ~isempty(config.paths.output)
    overallFigRoot    = fullfile(config.paths.output, 'figs',   'discprobe');
    overallTableRoot  = fullfile(config.paths.output, 'tables', 'discprobe');
else
    overallFigRoot    = fullfile(rootDir, 'output', 'figs',   'discprobe');
    overallTableRoot  = fullfile(rootDir, 'output', 'tables', 'discprobe');
end
if ~exist(overallFigRoot,   'dir'), mkdir(overallFigRoot);   end
if ~exist(overallTableRoot, 'dir'), mkdir(overallTableRoot); end

config.paths.globalDiscProbeFigRoot   = overallFigRoot;
config.paths.globalDiscProbeTableRoot = overallTableRoot;

doFlag = @(name,default) ...
    ( isfield(config,'do') && isfield(config.do,name) && ...
      all(logical(config.do.(name)(:))) ) ...
    || ( ~isfield(config,'do') && default );

for d = 1:numel(DATES)
    dateStr = char(DATES{d});
    fprintf('\n========== %s ===========\n', dateStr);
    
    % session-level paths
    SPaths = get_session_paths(config, dateStr);
    
    if ~isfolder(SPaths.session), mkdir(SPaths.session); end
    if ~isfolder(SPaths.units),   mkdir(SPaths.units);   end
    if ~isfolder(SPaths.tables),  mkdir(SPaths.tables);  end
    if ~isfolder(SPaths.figs),    mkdir(SPaths.figs);    end
    
    sessFigRoot    = fullfile(SPaths.figs,   'discprobe');
    sessTablesRoot = fullfile(SPaths.tables, 'discprobe');
    
    if ~isfolder(sessFigRoot),    mkdir(sessFigRoot);    end
    if ~isfolder(sessTablesRoot), mkdir(sessTablesRoot); end
    
    % Pass exact unit path to config so helpers can find it
    config.paths.figRoot    = sessFigRoot;
    config.paths.unitsRoot  = SPaths.units; 
    config.paths.tablesRoot = sessTablesRoot;
    
    % load session
    S = load_session(dateStr, config);

    % --- CRITICAL FIX ---
    % Force the Session ID/Date in the struct to match the FOLDER name.
    % This prevents U.dateStr from reverting to "250513" when the folder is "Jacomo_250513"
    S.dateStr   = dateStr; 
    S.sessionID = dateStr; 
    % --------------------
    
    % drop units with 0 spikes
    S = filter_zero_spike_units(S, config);
    
    % trial index table
    T = make_trial_index_table(S, COL, config);
    
    Tall  = table();
    Tpeak = table();
    Tcos  = table();
    
    for unitIdx = 1:S.nUnits
        % Safe ID printing
        if isfield(S, 'phyIDs'), currID = S.phyIDs(unitIdx); 
        elseif isfield(S, 'unitIDs'), currID = S.unitIDs(unitIdx); 
        else, currID = unitIdx; end

        fprintf('  unit %d/%d (ID=%d)\n', unitIdx, S.nUnits, currID);
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
            if isfield(stats, 'cosine') && ~isempty(stats.cosine)
                COS = stats.cosine;
                rowC = table;
                rowC.dateStr   = string(U.dateStr);
                rowC.unitID    = U.unitID;
                rowC.unitIdx   = unitIdx;
                rowC.unitType  = string(U.unitType);
                rowC.amp       = COS.amp;
                rowC.prefDeg   = COS.pref_deg;
                rowC.R2        = COS.R2;
                if isfield(COS, 'beta') && numel(COS.beta) >= 3
                    rowC.beta0   = COS.beta(1);
                    rowC.betaCos = COS.beta(2);
                    rowC.betaSin = COS.beta(3);
                else
                    rowC.beta0   = NaN;
                    rowC.betaCos = NaN;
                    rowC.betaSin = NaN;
                end
                Tcos = [Tcos; rowC]; %#ok<AGROW>
            end
        end
        if doFlag('anova', false)
            stats.anova = anova_analysis(U, config);
        end
        if doFlag('peakModel', false)
            stats.peak = peak_model(U, config);
        end
        
        % peak-model summary table
        if isfield(stats, 'peak') && ~isempty(stats.peak)
            P = stats.peak;
            t = table;
            t.dateStr        = string(U.dateStr);
            t.unitID         = U.unitID;
            t.unitIdx        = unitIdx;
            t.unitType       = string(U.unitType);
            if isfield(U, 'phyID')
                t.phyID      = U.phyID;
            else
                t.phyID      = U.unitID;
            end
            if isfield(P, 'class')
                t.class      = string(P.class);
            else
                t.class      = string(missing);
            end
            if isfield(P, 'amp1')
                t.amp1       = P.amp1;
            else
                t.amp1       = NaN;
            end
            if isfield(P, 'amp2')
                t.amp2       = P.amp2;
            else
                t.amp2       = NaN;
            end
            if isfield(P, 'amp_ratio')
                t.amp_ratio  = P.amp_ratio;
            else
                t.amp_ratio  = NaN;
            end
            if isfield(P, 'p_tuned')
                t.p_tuned    = P.p_tuned;
            else
                t.p_tuned    = NaN;
            end
            if isfield(P, 'p_bimodal')
                t.p_bimodal  = P.p_bimodal;
            else
                t.p_bimodal  = NaN;
            end
            if isfield(P, 'prefDeg')
                t.prefDeg    = P.prefDeg;
            elseif isfield(P, 'pref_deg')
                t.prefDeg    = P.pref_deg;
            elseif isfield(P, 'prefDirDeg')
                t.prefDeg    = P.prefDirDeg;
            else
                t.prefDeg    = NaN;
            end
            if isfield(P, 'FWHMDeg')
                t.FWHMDeg    = P.FWHMDeg;
            elseif isfield(P, 'FWHM_deg')
                t.FWHMDeg    = P.FWHM_deg;
            elseif isfield(P, 'fwhm_deg')
                t.FWHMDeg    = P.fwhm_deg;
            else
                t.FWHMDeg    = NaN;
            end
            Tpeak = [Tpeak; t]; %#ok<AGROW>
        end
        
        % === plotting ===
        if doFlag('rasterCanonical', true)
            plot_rasters(U, stats, COL, config, sessFigRoot, "canonical");
        end
        if doFlag('rasterByMean', false)
            plot_rasters(U, stats, COL, config, sessFigRoot, "mean");
        end
        if doFlag('tuningFigure', true)
            plot_tuning_curves(U, stats, COL, config, sessFigRoot);
        end
        if doFlag('rayleigh', false)
            plot_rayleigh(U, stats, COL, config, sessFigRoot);
        end
        if doFlag('fourierFig', false)
            plot_fourier(U, stats, COL, config, sessFigRoot);
        end
        if doFlag('boxPlots', false)
            plot_boxplots(U, stats, COL, config, sessFigRoot);
        end
        if doFlag('spikeAccounting', false)
            plot_spike_accounting(U, stats, COL, config, sessFigRoot);
        end
        if doFlag('dklSuite', true)
            plot_dkl_suite(U, stats, COL, config, sessFigRoot);
        end
        
        % per-unit .mat + quick CSV summary
        save_unit_outputs(config, U, stats);
        
        % simple per-unit summary row
        if isfield(stats, 'hueMeans') && isfield(stats.hueMeans, 'rate_mean')
            row = table;
            row.dateStr         = string(U.dateStr);
            row.unitID          = U.unitID;
            row.unitIdx         = unitIdx;
            row.unitType        = string(U.unitType);
            row.meanRateOverall = mean(stats.hueMeans.rate_mean, 'omitnan');
            Tall = [Tall; row]; %#ok<AGROW>
        end
    end
    
    % per-date unit summary CSV (per-session)
    if ~isempty(Tall)
        sumFile = fullfile(sessTablesRoot, sprintf('%s_unit_summary.csv', dateStr));
        try
            writetable(Tall, sumFile);
        catch ME
            warning('Could not write per-date summary CSV for %s: %s', dateStr, ME.message);
        end
        
        sumFileOverall = fullfile(overallTableRoot, sprintf('%s_unit_summary.csv', dateStr));
        try
            writetable(Tall, sumFileOverall);
        catch ME
            warning('Could not write overall unit summary CSV for %s: %s', dateStr, ME.message);
        end
    end
    
    % per-date peak summary CSV
    if ~isempty(Tpeak)
        peakFile = fullfile(sessTablesRoot, sprintf('%s_peak_summary.csv', dateStr));
        try
            writetable(Tpeak, peakFile);
        catch ME
            warning('Could not write peak summary CSV for %s: %s', dateStr, ME.message);
        end
        
        peakFileOverall = fullfile(overallTableRoot, sprintf('%s_peak_summary.csv', dateStr));
        try
            writetable(Tpeak, peakFileOverall);
        catch ME
            warning('Could not write overall peak summary CSV for %s: %s', dateStr, ME.message);
        end
    end
    
    % per-date cosine-fit summary CSV
    if ~isempty(Tcos)
        cosFile = fullfile(sessTablesRoot, sprintf('%s_cosine_summary.csv', dateStr));
        try
            writetable(Tcos, cosFile);
        catch ME
            warning('Could not write cosine summary CSV for %s: %s', dateStr, ME.message);
        end
        
        cosFileOverall = fullfile(overallTableRoot, sprintf('%s_cosine_summary.csv', dateStr));
        try
            writetable(Tcos, cosFileOverall);
        catch ME
            warning('Could not write overall cosine summary CSV for %s: %s', dateStr, ME.message);
        end
    end
end
fprintf('\nDone.\n');
end