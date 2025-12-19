function run_vss_peakmodels(DATES)
% RUN_VSS_PEAKMODELS - Wrapper for VSS abstract analysis.
% Fits harmonic peak models (unimodal vs bimodal) to max saturation trials.
% Generates CSV tables for peak classification and hue means.
% SAVES TO: 
%   1. output/vss/<Folder_Name>/tables/
%   2. data/<Folder_Name>/tables/vss/ (Local copy)

    clc;
    % --- 1. SETUP & CONFIG ---
    rootDir = fileparts(mfilename('fullpath'));
    addpath(genpath(rootDir));
    
    config = load_config(); 
    COL    = load_color_metadata(config); 
    
    % --- 2. DATE HANDLING ---
    if nargin < 1 || isempty(DATES)
        DATES = 'all';
    end
    
    % Expand 'all' into list of session folders
    if ischar(DATES) || isstring(DATES)
        ds = strtrim(char(DATES));
        if strcmpi(ds, 'all')
            baseRoot = config.paths.base; 
            dInfo = dir(baseRoot);
            names = {dInfo.name};
            % Find folders that look like dates (6 digits) or SessionIDs (Name_Date)
            isDateDir = [dInfo.isdir] & ~ismember(names, {'.','..'});
            % Filter for folders containing the date pattern
            validIdx = find(isDateDir);
            finalList = {};
            for k = validIdx
                if regexp(names{k}, '\d{6}')
                    finalList{end+1} = names{k}; %#ok<AGROW>
                end
            end
            DATES = sort(finalList);
        else
            DATES = {ds};
        end
    end
    
    % --- 3. FLAGS & SEEDS ---
    if isfield(config, 'pt') && isfield(config.pt, 'seed')
        rng(config.pt.seed);
    end
    
    doFlag = @(name,def) (isfield(config,'analysis') && isfield(config.analysis,name) && ...
                         all(logical(config.analysis.(name)(:)))) || ...
                         (~isfield(config,'analysis') && def);
                     
    DO_DROP_OUTLIERS   = doFlag('dropOutliers', true);
    DO_SAVE_UNIT_MATS  = doFlag('saveVssUnitMats', true);
    DO_BUILD_KEEP_LIST = doFlag('buildVssKeep', true);
    
    allKeepRows = table();

    % --- 4. OUTPUT PATHS (Centralized) ---
    vssRoot = fullfile(config.paths.output, 'vss');
    if ~exist(vssRoot, 'dir'), mkdir(vssRoot); end
    
    globalTableDir = fullfile(config.paths.output, 'tables');
    if ~exist(globalTableDir, 'dir'), mkdir(globalTableDir); end

    % --- 5. SESSION LOOP ---
    for d = 1:numel(DATES)
        % DATES{d} might be '250513' or 'Jacomo_250513'
        % We need to resolve the exact folder name on disk
        targetDateStr = char(DATES{d});
        
        % Helper to find the actual folder name (e.g. transforms '250513' -> 'Jacomo_250513')
        sessionFolderName = resolve_session_folder(config.paths.base, targetDateStr);
        
        if isempty(sessionFolderName)
            warning('Folder for %s not found in %s. Skipping.', targetDateStr, config.paths.base);
            continue;
        end
        
        fprintf('\n=== VSS Peak Model Analysis: %s ===\n', sessionFolderName);
        
        % ---------------- PATH SETUP ----------------
        % 1. Central Output Path: output/vss/<SessionFolder>/tables
        vssSessDir = fullfile(vssRoot, sessionFolderName);
        vssTables  = fullfile(vssSessDir, 'tables');
        vssMats    = fullfile(vssSessDir, 'units');
        
        if ~exist(vssSessDir, 'dir'), mkdir(vssSessDir); end
        if ~exist(vssTables, 'dir'),  mkdir(vssTables);  end
        if DO_SAVE_UNIT_MATS && ~exist(vssMats, 'dir'), mkdir(vssMats); end
        
        % 2. Local Session Path: data/<SessionFolder>/tables/vss
        localSessDir   = fullfile(config.paths.base, sessionFolderName);
        localTablesDir = fullfile(localSessDir, 'tables', 'vss');
        if ~exist(localTablesDir, 'dir'), mkdir(localTablesDir); end
        % --------------------------------------------

        % Load Session (Use the folder name or date string)
        try
            S = load_session(sessionFolderName, config);
        catch ME
            warning('Could not load session %s: %s', sessionFolderName, ME.message);
            continue;
        end
        
        S = filter_zero_spike_units(S, config);
        T = make_trial_index_table(S, COL, config); 
        
        peakRows = table();
        meanRows = table();
        
        alpha = 0.05;
        if isfield(config,'tuning') && isfield(config.tuning,'alphaHarmonic')
            alpha = config.tuning.alphaHarmonic;
        end
        
        % --- UNIT LOOP ---
        for unitIdx = 1:S.nUnits
            if unitIdx > numel(S.unitIDs), break; end
            uID = S.unitIDs(unitIdx);
            
            if isfield(S, 'unitTypes') && numel(S.unitTypes) >= unitIdx
                uType = S.unitTypes{unitIdx};
            else
                uType = '?';
            end
            
            fprintf('  Processing Unit %d/%d (ID=%d %s)...\n', unitIdx, S.nUnits, uID, uType);
            
            U = struct();
            U.dateStr  = sessionFolderName; % Ensure we use the full folder name
            U.unitID   = uID;
            U.unitType = uType;
            U.spk      = S.spk{unitIdx}; 
            U.trials   = T;              
            U.nTrials  = height(T);
            
            if isfield(config,'window') && isfield(config.window,'early')
                U.winEarly = config.window.early;
            else
                U.winEarly = [0.05 0.25]; 
            end

            if DO_DROP_OUTLIERS
                U = detect_outliers(U, config);
            end
            
            U_high = subset_unit_to_max_saturation(U, 1.0);
            
            if ~isfield(U_high,'nTrials') || U_high.nTrials < 10
                continue; 
            end
            
            hm = compute_hue_means(U_high);      
            P  = peak_model_harmonic(U_high, config); 
            
            dateStrStr   = string(sessionFolderName);
            unitTypeStr  = string(U.unitType);
            unitIDVal    = U.unitID;
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
            
            nH = numel(hm.hues);
            if nH > 0
                newMean = table( ...
                    repmat(dateStrStr,nH,1), repmat(unitTypeStr,nH,1), repmat(unitIDVal,nH,1), ...
                    repmat(unitLabelStr,nH,1), hm.hues(:), hm.rate_mean(:), hm.rate_sem(:), hm.n_trials(:), ...
                    'VariableNames', {'dateStr','unitType','unitID','unitLabel','hueID','rate_mean','rate_sem','n_trials'});
                meanRows = [meanRows; newMean];
            end
            
            peakRow = table(dateStrStr, unitTypeStr, unitIDVal, unitLabelStr, ...
                string(P.class), string(P.class), isTunedFit, kMax_fit, fracMax_fit, ...
                P.amp1, P.amp2, P.amp_ratio, nH, U_high.nTrials, ...
                'VariableNames',{'dateStr','unitType','unitID','unitLabel','peakClass','peakClass_orig',...
                'isTunedFit','kMax_fit','fracMax_fit','amp1','amp2','amp_ratio','nHues','nTrials'});
            
            if isfield(P,'R2'), peakRow.R2 = P.R2; else, peakRow.R2 = NaN; end
            if isfield(P,'p_tuned'), peakRow.p_tuned = P.p_tuned; else, peakRow.p_tuned = NaN; end
            if isfield(P,'p_bimodal'), peakRow.p_bimodal = P.p_bimodal; else, peakRow.p_bimodal = NaN; end
            
            peakRows = [peakRows; peakRow];
            
            if DO_SAVE_UNIT_MATS
                saveName = fullfile(vssMats, sprintf('%s_%s%d_VSS_Peak.mat', sessionFolderName, U.unitType, U.unitID));
                save(saveName, 'U_high', 'hm', 'P');
            end
            
            isKeep = ismember(P.class, ["unimodal","bimodal"]);
            if DO_BUILD_KEEP_LIST && isKeep
                allKeepRows = [allKeepRows; table( ...
                    dateStrStr, unitLabelStr, 1, string(P.class), unitTypeStr, unitIDVal, ...
                    'VariableNames', {'session','unit','keep','peakClass','unitType','unitID'})];
            end
        end
        
        % --- SAVE SESSION TABLES ---
        meanName = sprintf('%s_vss_hueMeans_sat1.csv', sessionFolderName);
        peakName = sprintf('%s_vss_peakModel_sat1.csv', sessionFolderName);
        
        if ~isempty(meanRows)
            writetable(meanRows, fullfile(vssTables, meanName));
            writetable(meanRows, fullfile(localTablesDir, meanName)); 
        end
        if ~isempty(peakRows)
            writetable(peakRows, fullfile(vssTables, peakName));
            writetable(peakRows, fullfile(localTablesDir, peakName)); 
        end
        
        fprintf('  > Saved session tables to: %s\n', vssTables);
        fprintf('  > Saved local copies to:   %s\n', localTablesDir);
    end
    
    % --- 6. SAVE GLOBAL KEEP LIST ---
    if DO_BUILD_KEEP_LIST && ~isempty(allKeepRows)
        keepFile = fullfile(globalTableDir, 'ALL_keep_vss.csv');
        writetable(allKeepRows, keepFile);
        fprintf('\n✅ Wrote Global VSS Keep List to: %s\n', keepFile);
        fprintf('   Total units kept: %d\n', height(allKeepRows));
    end
end

% -------------------------------------------------------------------------
% local helpers
% -------------------------------------------------------------------------
function U_out = subset_unit_to_max_saturation(U, maxSatVal)
    U_out = U;
    if ~isfield(U,'trials'), U_out.nTrials=0; return; end
    
    if ismember('satID', U.trials.Properties.VariableNames)
        satCol = U.trials.satID;
    elseif ismember('satVal', U.trials.Properties.VariableNames)
        satCol = U.trials.satVal;
    elseif ismember('sat', U.trials.Properties.VariableNames)
        satCol = U.trials.sat;
    else
        warning('Cannot find saturation column (satID/satVal) in U.trials');
        U_out.nTrials = 0;
        return;
    end
    
    mask = abs(satCol - maxSatVal) < 0.01;
    
    U_out.trials = U.trials(mask, :);
    if isfield(U,'spk') && numel(U.spk) == U.nTrials
        U_out.spk = U.spk(mask);
    end
    U_out.nTrials = sum(mask);
end