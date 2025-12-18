function run_wachtler_analysis(DATES, MODE)
% Wachtler-style color tuning analysis using main DiscProbe pipeline helpers.

clc;

if nargin < 1
    DATES = [];
end
if nargin < 2 || isempty(MODE)
    MODE = 'full';  % 'full' or 'fitsOnly'
end


% make sure repo is on path
rootDir = fileparts(mfilename('fullpath'));
addpath(genpath(rootDir));

% core config + color metadata
config = load_config();
COL    = load_color_metadata(config);
COL    = add_lms_axes(COL);   % uses COL.lms to build L-M / S / Lum

monkey = config.monkey;

% handle DATES / "all" 
if nargin < 1
    DATES = [];
end

if ischar(DATES) || isstring(DATES)
    if strcmpi(strtrim(DATES), 'all')
        DATES = [];
    end
end

PATHS = struct();
PATHS.baseDiscProbeLocal = config.paths.base;
PATHS.baseDiscProbesCode = config.paths.code;
PATHS.baseOut            = config.paths.output;

if isempty(DATES)
    pat = sprintf('%s_*_ExpTrialsDisc.mat', monkey);
    D = dir(fullfile(PATHS.baseDiscProbeLocal, '**', pat));
    if isempty(D)
        error('run_wachtler_analysis:NoSessions', ...
            'No ExpTrialsDisc files found under %s', PATHS.baseDiscProbeLocal);
    end

    dateList = strings(0,1);
    for k = 1:numel(D)
        % look at the *folder* that contains the file
        [~, parentName] = fileparts(D(k).folder);

        % only keep if the folder name is exactly 6 digits, e.g. '250513'
        tok = regexp(parentName, '^\d{6}$', 'match', 'once');
        if ~isempty(tok)
            dateList(end+1,1) = string(tok);
        end
    end

    dateList = unique(dateList);
    if isempty(dateList)
        error('run_wachtler_analysis:NoDateFolders', ...
            'Found ExpTrialsDisc files, but none in 6-digit date folders under %s', ...
            PATHS.baseDiscProbeLocal);
    end
else
    dateList = normalize_dates(DATES);
    dateList = string(dateList(:));
end

% ---------- Wachtler-specific settings ----------
config.tuning.win         = [0.05 0.15];
config.tuning.baselineWin = [-0.15 0];
config.tuning.useSaturIDs = [0.2 0.33 0.5 0.66 0.8 1];

if ~isfield(config, 'tuning') || ~isfield(config.tuning, 'R2_thresh')
    config.tuning.R2_thresh = 0.5;
end

isFitsOnly = strcmpi(MODE, 'fitsOnly');

if isFitsOnly
    config.doPlots.tuningPerUnit    = false;
    config.doPlots.tuningPopulation = false;
else
    config.doPlots.tuningPerUnit    = true;
    config.doPlots.tuningPopulation = true;
end


% ---------- output roots ----------
if isfield(config, 'paths') && isfield(config.paths, 'output') && ~isempty(config.paths.output)
    outRoot = config.paths.output;
else
    outRoot = fullfile(rootDir, 'output');
end

wachtRoot = fullfile(outRoot, 'wachtler');
if ~isfolder(wachtRoot), mkdir(wachtRoot); end

for d = 1:numel(dateList)
    dateStr = char(dateList(d));
    fprintf('\n========== Wachtler analysis %s ==========\n', dateStr);

    sessionRoot = fullfile(wachtRoot, dateStr);
    figRoot     = fullfile(sessionRoot, 'figs');
    tablesRoot  = fullfile(sessionRoot, 'tables');

    if ~isfolder(sessionRoot), mkdir(sessionRoot); end
    if ~isfolder(figRoot),     mkdir(figRoot);     end
    if ~isfolder(tablesRoot),  mkdir(tablesRoot);  end

    config.paths.wachtlerFigRoot    = figRoot;
    config.paths.wachtlerTablesRoot = tablesRoot;

    % ---------- load and prep session ----------
    S = load_session(dateStr, config);
    S = filter_zero_spike_units(S, config);

    T = make_trial_index_table(S, COL, config);

    nUnits = S.nUnits;
    allUnitTuning = cell(nUnits, 1);

    for unitIdx = 1:nUnits
        fprintf('  unit %d/%d (ID=%d)\n', unitIdx, nUnits, S.unitIDs(unitIdx));

        U = prepare_unit(S, unitIdx, T, config);
        U = attach_color_index(U, COL);  % maps trials to color rows

        if isfield(config,'do') && isfield(config.do,'dropOutliers') && config.do.dropOutliers
            U = detect_outliers(U, config);
        end

        unitTuning = compute_unit_color_tuning_wachtler(U, COL, config);
        allUnitTuning{unitIdx} = unitTuning;

        doPerUnit = isfield(config, 'doPlots') && ...
                    isfield(config.doPlots, 'tuningPerUnit') && ...
                    config.doPlots.tuningPerUnit;

        if doPerUnit && ~isempty(unitTuning)
            % classic Wachtler per-unit figure
            H = plot_unit_color_tuning_wachtler(U, unitTuning, COL, config);
            if isfield(H,'fig') && ishghandle(H.fig)
                tag = sprintf('%s_%s_%d', string(U.dateStr), U.unitType, U.unitID);
                outPng = fullfile(figRoot, [tag '_wachtler_tuning.png']);
                try
                    saveas(H.fig, outPng);
                catch ME
                    warning('Could not save Wachtler tuning fig for unit %d: %s', ...
                        U.unitID, ME.message);
                end
                close(H.fig);
            end

            % new: 3D LMS blob for this unit
            tuning3D = compute_unit_tuning_lms(U, COL, config);
            if ~isempty(tuning3D)
                H3 = plot_unit_3d_blob_lms(U, tuning3D, config);
                if isfield(H3,'fig') && ishghandle(H3.fig)
                    tag = sprintf('%s_%s_%d', string(U.dateStr), U.unitType, U.unitID);
                    outPng = fullfile(figRoot, [tag '_wachtler_blob3D.png']);
                    try
                        saveas(H3.fig, outPng);
                    catch ME
                        warning('Could not save Wachtler 3D blob for unit %d: %s', ...
                            U.unitID, ME.message);
                    end
                    close(H3.fig);
                end
            end
        end
    end

    % ---------- build population summary + CSV with fit stats ----------
    if ~isempty(allUnitTuning)
        unitID  = S.unitIDs(:);
        unitNum = (1:numel(unitID)).';

        nU = numel(allUnitTuning);
        R2      = nan(nU,1);
        chi2val = nan(nU,1);
        df      = nan(nU,1);
        pChi2   = nan(nU,1);
        prefDeg = nan(nU,1);
        fwhmDeg = nan(nU,1);

        for ii = 1:nU
            ut = allUnitTuning{ii};
            if isempty(ut) || ~isfield(ut, 'fit') || isempty(ut.fit)
                continue;
            end
            f = ut.fit;

            if isfield(f,'R2'),        R2(ii)      = f.R2;        end
            if isfield(f,'chi2'),      chi2val(ii) = f.chi2;      end
            if isfield(f,'df'),        df(ii)      = f.df;        end
            if isfield(f,'pChi2'),     pChi2(ii)   = f.pChi2;     end
            if isfield(f,'phi0_deg'),  prefDeg(ii) = f.phi0_deg;  end
            if isfield(f,'fwhm_deg'),  fwhmDeg(ii) = f.fwhm_deg;  end
        end

        r2Thresh  = config.tuning.R2_thresh;
        isTunedR2 = R2 >= r2Thresh;

        nMin    = min([numel(unitID), numel(R2)]);
        unitID  = unitID(1:nMin);
        unitNum = unitNum(1:nMin);
        R2      = R2(1:nMin);
        chi2val = chi2val(1:nMin);
        df      = df(1:nMin);
        pChi2   = pChi2(1:nMin);
        prefDeg = prefDeg(1:nMin);
        fwhmDeg = fwhmDeg(1:nMin);
        isTunedR2 = isTunedR2(1:nMin);

        Tfits = table(unitID, unitNum, R2, chi2val, df, pChi2, ...
                      isTunedR2, prefDeg, fwhmDeg, ...
            'VariableNames', {'unitID','unitNum','R2','chi2','df','pChi2', ...
                              'isTunedR2','prefDeg','FWHMDeg'});

        csvFile = fullfile(tablesRoot, sprintf('%s_wachtler_fits.csv', dateStr));
        try
            writetable(Tfits, csvFile);
            fprintf('  Saved Wachtler fits CSV: %s\n', csvFile);
        catch ME
            warning('Could not write Wachtler fits CSV for %s: %s', dateStr, ME.message);
        end

        doPop = isfield(config, 'doPlots') && ...
                isfield(config.doPlots, 'tuningPopulation') && ...
                config.doPlots.tuningPopulation;

        if doPop
            pop = summarize_tuning_population(allUnitTuning, config);
            Hpop = plot_tuning_population_histograms(pop, config);
            if isfield(Hpop,'fig') && ishghandle(Hpop.fig)
                outPng = fullfile(sessionRoot, sprintf('%s_wachtler_population.png', dateStr));
                try
                    saveas(Hpop.fig, outPng);
                catch ME
                    warning('Could not save Wachtler population fig for %s: %s', ...
                        dateStr, ME.message);
                end
                close(Hpop.fig);
            end
        end
    end
end

fprintf('\nDone.\n');

end
