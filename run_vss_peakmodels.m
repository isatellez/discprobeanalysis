function run_vss_peakmodels(DATES)

%this script is the wrapper i used for the data for abstract submission for vss

clc;


%the file lives in the root directory
rootDir = fileparts(mfilename('fullpath'));
addpath(genpath(rootDir));

config = load_config(); %load config
COL    = load_color_metadata(config); %prepare color look up tables etc

if nargin < 1 || isempty(DATES)
    DATES = {'all'};
end
%by default, the wrapper will run every session

% for 'all' option 
if ischar(DATES) || isstring(DATES)
    ds = strtrim(char(DATES));
    if strcmpi(ds, 'all')
        baseRoot = config.paths.base; %where mats live
        dInfo = dir(baseRoot);
        names = {dInfo.name};
        isDateDir = [dInfo.isdir] & ...
            ~ismember(names, {'.','..'}) & ...
            ~cellfun(@isempty, regexp(names, '^\d{6}$', 'once'));
        DATES = sort(names(isDateDir));
    else
        DATES = {ds};
    end
end

DATES = normalize_dates(DATES);

% RNG & flags
if isfield(config, 'pt') && isfield(config.pt, 'seed')
    rng(config.pt.seed);
end

doFlag = @(name,def) (isfield(config,'do') && isfield(config.do,name) && ...
                     all(logical(config.do.(name)(:)))) || ...
                     (~isfield(config,'do') && def);

%do you want to drop outliers? (usually ALWAYS yes)
DO_DROP_OUTLIERS   = doFlag('dropOutliers',true);
%do you want to save mats?
DO_SAVE_UNIT_MATS  = doFlag('saveVssUnitMats',true);
%do you want to build the "keep" list that is used in other analysis?
DO_BUILD_KEEP_LIST = doFlag('buildVssKeep',true);

%empty table that will have all out kept units
allKeepRows = table();

for d = 1:numel(DATES)
    dateStr = char(DATES{d});
    fprintf('\n single vs two peak-model analysis: %s\n', dateStr);

    %get session paths so it knows where to save stuff
    SPaths = get_session_paths(config, dateStr);
    if ~exist(SPaths.session,'dir'), mkdir(SPaths.session); end
    if ~exist(SPaths.units,'dir'),   mkdir(SPaths.units);   end
    if ~exist(SPaths.tables,'dir'),  mkdir(SPaths.tables);  end

    % lets store VSS results under /vss/<date>/tables
    vssRoot = fullfile(config.paths.output, 'vss');
    vssSess = fullfile(vssRoot, dateStr); 
    vssTables = fullfile(vssSess, 'tables');
    if ~exist(vssTables, 'dir'), mkdir(vssTables); end

    %load .mat for that session with all the spikes 
    try
        S = load_session(dateStr, config);
    catch ME
        warning('Could not load session %s: %s', dateStr, ME.message);
        continue;
    end

    %preparing unit
    S = filter_zero_spike_units(S, config); %drop units with 0 spikes
    T = make_trial_index_table(S, COL, config); %make trial index with all trials for the session
    %this is to get their rgbs, dkl, lms values per trial

    peakRows = table();
    meanRows = table();
    fourRows = table();

    alpha = 0.05;
    if isfield(config,'tuning') && isfield(config.tuning,'alphaHarmonic')
        alpha = config.tuning.alphaHarmonic;
    end

    for unitIdx = 1:S.nUnits
        fprintf('  unit %d/%d (ID=%d)\n', unitIdx, S.nUnits, S.unitIDs(unitIdx));
        U = prepare_unit(S, unitIdx, T, config);
        %get all spike info for unit and unit metadata


        %drop outliers (spikes that fall beyond 1.5 IQR for every trial
        %type
        if DO_DROP_OUTLIERS
            U = detect_outliers(U, config);
        end


        % only max saturation trials
        U_high = subset_unit_to_max_saturation(U, 1.0);
        if ~isfield(U_high,'nTrials') || U_high.nTrials == 0
            continue;
        end

        hm = compute_hue_means(U_high);
        P  = peak_model_harmonic(U_high, config);   % <-- fit-based classifier
        FT = fourier_analysis(U_high, config);      % descriptive only

        dateStrStr   = string(U_high.dateStr);
        unitTypeStr  = string(U_high.unitType);
        unitIDVal    = U_high.unitID;
        unitLabelStr = unitTypeStr + string(unitIDVal);

        isTunedFit = isfield(P,'p_tuned') && ~isnan(P.p_tuned) && (P.p_tuned < alpha);

        if P.class == "unimodal"
            kMax_fit = 1;
        elseif P.class == "bimodal"
            kMax_fit = 2;
        else
            kMax_fit = NaN;
        end

        fracMax_fit = NaN;

        % -------- hue mean table --------
        nH = numel(hm.hues);
        if nH > 0
            meanRows = [meanRows; table( ...
                repmat(dateStrStr,nH,1), repmat(unitTypeStr,nH,1), repmat(unitIDVal,nH,1), ...
                repmat(unitLabelStr,nH,1), hm.hues(:), hm.rate_mean(:), hm.rate_sem(:), hm.n_trials(:), ...
                'VariableNames', {'dateStr','unitType','unitID','unitLabel','hueID','rate_mean','rate_sem','n_trials'})];
        end

        % -------- peak model table --------
        peakRow = table(dateStrStr,unitTypeStr,unitIDVal,unitLabelStr,...
            string(P.class),string(P.class),isTunedFit,kMax_fit,fracMax_fit,...
            P.amp1,P.amp2,P.amp_ratio,numel(hm.hues),U_high.nTrials,...
            'VariableNames',{'dateStr','unitType','unitID','unitLabel','peakClass','peakClass_orig',...
            'isTunedFit','kMax_fit','fracMax_fit','amp1','amp2','amp_ratio','nHues','nTrials'});

       if isfield(P,'R2')
        peakRow.R2 = P.R2;
        else
            peakRow.R2 = NaN;
        end
        
        if isfield(P,'p_tuned')
            peakRow.p_tuned = P.p_tuned;
        else
            peakRow.p_tuned = NaN;
        end
        
        if isfield(P,'p_bimodal')
            peakRow.p_bimodal = P.p_bimodal;
        else
            peakRow.p_bimodal = NaN;
        end

        peakRows = [peakRows; peakRow];

        if DO_SAVE_UNIT_MATS
            save_vss_unit(config, U_high, hm, P);
        end

     % --- decide if this unit qualifies for the global keep list ---
isKeep = ismember(P.class, ["unimodal","bimodal"]);

if DO_BUILD_KEEP_LIST && isKeep
    allKeepRows = [allKeepRows; table( ...
        dateStrStr, unitLabelStr, 1, string(P.class), unitTypeStr, unitIDVal, ...
        'VariableNames', {'session','unit','keep','peakClass','unitType','unitID'})];
end


    end

    % --- write per-session tables under /vss/<date>/tables ---
    if ~isempty(meanRows)
        writetable(meanRows, fullfile(vssTables, sprintf('%s_vss_hueMeans_sat1.csv', dateStr)));
    end
    if ~isempty(peakRows)
        writetable(peakRows, fullfile(vssTables, sprintf('%s_vss_peakModel_sat1.csv', dateStr)));
    end
    if ~isempty(fourRows)
        writetable(fourRows, fullfile(vssTables, sprintf('%s_vss_fourier_sat1.csv', dateStr)));
    end
end

% --- global keeper ---
if DO_BUILD_KEEP_LIST && ~isempty(allKeepRows)
    keepFile = fullfile(config.paths.base, 'ALL_keep_vss.csv');
    writetable(allKeepRows, keepFile);
    fprintf('✅ Wrote global keep list to %s (%d units kept)\n', keepFile, height(allKeepRows));
end

end
