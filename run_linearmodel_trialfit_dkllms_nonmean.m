function run_linearmodel_trialfit_dkllms_nonmean(DATES)
% DiscProbe_STA_LMSLINv6_TrialFit - Fits LMS/DKL models to raw trial data.
% Fits models using Grand STA weights projected onto trial stimuli.
% This version uses trial-by-trial variance instead of stimulus means.
%
% SAVES TO:
%   1. output/sta/<Session>/ (CSV summary)
%   2. data/<Session>/sta/   (Local copy)
%   3. output/tables/        (Global summary)
    clc;
    rootDir = fileparts(mfilename('fullpath'));
    addpath(genpath(rootDir));
    
    config = load_config();
    COL    = load_color_metadata(config); 
    
    % --- 1. DATE HANDLING ---
    if nargin < 1 || isempty(DATES)
        DATES = 'all';
    end
    
    if ischar(DATES) || isstring(DATES)
        ds = strtrim(char(DATES));
        if strcmpi(ds, 'all')
            d = dir(config.paths.base);
            names = {d([d.isdir]).name};
            % Filter for folders looking like dates (6 digits) or SessionIDs
            validIdx = cellfun(@(s) ~isempty(regexp(s, '\d{6}', 'once')), names) & ...
                       ~ismember(names, {'.','..'});
            DATES = sort(names(validIdx));
        else
            DATES = {ds};
        end
    else
        DATES = cellstr(DATES);
    end
    warning('off','stats:glmfit:IllConditioned');
    warning('off','stats:glmfit:IterationLimit');
    warning('off','MATLAB:rankDeficientMatrix');
    timeWin     = [0 0.2];
    % Currently looking at all standard sat values: 0.33, 0.66, and 1.0
    saturIDs    = [0.33 0.66 1];
    nSat        = numel(saturIDs);
    N_BOOT_REPS = 1000;
    
    allRowsGlobal = [];
    % ------------- Keep List ----------------
    % Pass the OUTPUT/TABLES directory explicitly, or rely on the robust search in the function
    % We try to pass the exact location of the global tables first.
    keepMapDir = fullfile(config.paths.output, 'tables'); 
    keepMap = load_keep_map(keepMapDir); 
    
    if isempty(keepMap)
        % Fallback to base if empty (legacy support)
        keepMap = load_keep_map(config.paths.base);
    end
    % ------------- Color Metadata ----------------
    lmsFile = fullfile(config.paths.code, 'lms.csv');
    if ~isfile(lmsFile), error('Missing lms.csv at %s', lmsFile); end
    LMS_all = readmatrix(lmsFile);
    
    if isfield(COL, 'nHue')
        nHue = COL.nHue;
    else
        colsFile = fullfile(config.paths.code, 'DiscProbeColsUpdated.mat');
        if isfile(colsFile)
            tmp = load(colsFile, 'ProbeColIDs');
            nHue = max(tmp.ProbeColIDs(:,1));
        else
            nHue = 16; 
        end
    end
    
    try
        LMS_satur2 = reshape(LMS_all(1:nSat*nHue,1:3), [nSat, nHue, 3]);
    catch
        error('lms.csv rows (%d) do not match nSat*nHue (%d*%d)', size(LMS_all,1), nSat, nHue);
    end
    % --- OUTPUT PATHS ---
    outRootCentral = fullfile(config.paths.output, 'sta');
    if ~exist(outRootCentral, 'dir'), mkdir(outRootCentral); end
    
    outRootGlobal = fullfile(config.paths.output, 'tables');
    if ~exist(outRootGlobal, 'dir'), mkdir(outRootGlobal); end
    % =====================================================================
    % LOOP OVER SESSIONS
    % =====================================================================
    for idate = 1:numel(DATES)
        targetDateStr = char(DATES{idate});
        
        % Resolve actual folder name (e.g., '250513' -> 'Jacomo_250513')
        sessionFolderName = resolve_session_folder(config.paths.base, targetDateStr);
        if isempty(sessionFolderName)
            warning('Folder for %s not found. Skipping.', targetDateStr);
            continue;
        end
        
        fprintf('\n=== %s: Trial-by-Trial LMS & DKL fitting ===\n', sessionFolderName);
        % --- PATH SETUP ---
        localSessDir  = fullfile(config.paths.base, sessionFolderName);
        
        % 1. Central Output
        outDirCentral = fullfile(outRootCentral, sessionFolderName);
        if ~exist(outDirCentral, 'dir'), mkdir(outDirCentral); end
        
        % 2. Local Output
        outDirLocal = fullfile(localSessDir, 'sta');
        if ~exist(outDirLocal, 'dir'), mkdir(outDirLocal); end
        % --- LOAD SESSION ---
        try
            S = load_session(sessionFolderName, config);
        catch ME
            fprintf('  [%s] Load failed: %s\n', sessionFolderName, ME.message);
            continue;
        end
        
        if isfield(S, 'ExptTrialsDisc')
             T = make_trial_index_table(S, COL, config);
        else
             warning('ExptTrialsDisc missing in S.');
             continue;
        end
        
        nTrials = height(T);
        
        % Robust Column Handling
        if ismember('hueID', T.Properties.VariableNames), hueIdx = T.hueID;
        elseif ismember('hueIdx', T.Properties.VariableNames), hueIdx = T.hueIdx;
        elseif ismember('hue', T.Properties.VariableNames), hueIdx = T.hue;
        else, warning('T table missing hueID. Skipping.'); continue; end
        
        if ismember('satID', T.Properties.VariableNames), satVal = T.satID;
        elseif ismember('satVal', T.Properties.VariableNames), satVal = T.satVal;
        elseif ismember('sat', T.Properties.VariableNames), satVal = T.sat;
        else, satVal = nan(nTrials,1); end
        
        [satIdx, okSat] = ismember(satVal, saturIDs);
        validStim = okSat & hueIdx >= 1 & hueIdx <= nHue;
        
        % LMS Matrix
        trial_LMS = nan(nTrials, 3);
        for t = 1:nTrials
            if validStim(t)
                trial_LMS(t,:) = squeeze(LMS_satur2(satIdx(t), hueIdx(t), :));
            end
        end
        % Count Spikes
        spkcounts_all = nan(nTrials, S.nUnits);
        for cc = 1:S.nUnits
            unitSpikes = S.spk{cc};
            for t = 1:nTrials
                spks = unitSpikes{t};
                n = sum(spks >= timeWin(1) & spks < timeWin(2));
                spkcounts_all(t, cc) = n;
            end
        end
        allRows = [];
        
        % --- LOOP OVER UNITS ---
        for cc = 1:S.nUnits
            phyID = S.unitIDs(cc);
            if isfield(S, 'unitTypes')
                unitType = S.unitTypes{cc};
            else
                unitType = '?'; 
            end
            unitLabel = sprintf('%s%d', unitType, phyID);
            % --- CHECK KEEP LIST ---
            if ~isempty(keepMap)
                keyFull = sprintf('%s_%s', sessionFolderName, unitLabel);
                keyDate = sprintf('%s_%s', targetDateStr, unitLabel);
                
                keepVal = 0; 
                if isKey(keepMap, keyFull)
                    keepVal = keepMap(keyFull);
                elseif isKey(keepMap, keyDate)
                    keepVal = keepMap(keyDate);
                else
                    keepVal = 0; 
                end
                
                if keepVal == 0
                    continue; 
                end
            end
            
            % Locate STA
            staFileInfo = find_sta_file(localSessDir, unitType, phyID);
            if isempty(staFileInfo), continue; end
            
            try
                [w_chan, chanNames_sta] = grand_sta_per_color(staFileInfo.fullpath);
            catch ME
                fprintf('  [%s] %s: STA load failed: %s\n', sessionFolderName, unitLabel, ME.message);
                continue;
            end
            
            idxL = find(strcmpi(chanNames_sta,'L'),1);
            idxM = find(strcmpi(chanNames_sta,'M'),1);
            idxS = find(strcmpi(chanNames_sta,'S'),1);
            
            if isempty(idxL) || isempty(idxM), continue; end
            
            wL = w_chan(idxL); 
            wM = w_chan(idxM);
            hasS = ~isempty(idxS); 
            if hasS, wS = w_chan(idxS); end
            
            L_disc_all = trial_LMS(:,1);
            M_disc_all = trial_LMS(:,2);
            if hasS, S_disc_all = trial_LMS(:,3); end
            
            y_all = spkcounts_all(:, cc);
            
            if hasS
                valid = validStim & isfinite(L_disc_all) & isfinite(M_disc_all) & isfinite(S_disc_all) & isfinite(y_all);
            else
                valid = validStim & isfinite(L_disc_all) & isfinite(M_disc_all) & isfinite(y_all);
            end
            
            if nnz(valid) < 50, continue; end
            
            y = y_all(valid);
            L_disc = L_disc_all(valid);
            M_disc = M_disc_all(valid);
            if hasS, S_disc = S_disc_all(valid); end
            
            % --- LMS Trial Fit ---
            gL = wL * L_disc; gM = wM * M_disc;
            if hasS
                gS = wS * S_disc;
                G_LMS = [gL, gM, gS]; chanNames_LMS = {'L','M','S'};
            else
                G_LMS = [gL, gM]; chanNames_LMS = {'L','M'};
            end
            
            [~, bestChan_LMS, LMS_b0, LMS_betas, R2_best_trial_LMS, R2_all_trial_LMS] = ...
                compute_space_stats_trial(G_LMS, y, chanNames_LMS);
            
            LMS_bL=NaN; LMS_bM=NaN; LMS_bS=NaN;
            for ii=1:numel(chanNames_LMS)
                switch upper(chanNames_LMS{ii})
                    case 'L', LMS_bL=LMS_betas(ii);
                    case 'M', LMS_bM=LMS_betas(ii);
                    case 'S', LMS_bS=LMS_betas(ii);
                end
            end
            
            % --- DKL Trial Fit ---
            bestChan_DKL=''; R2_best_trial_DKL=NaN; R2_all_trial_DKL=NaN;
            DKL_b0=NaN; DKL_bLmM=NaN; DKL_bSLpM=NaN;
            
            if hasS
                g_LmM = wL*L_disc - wM*M_disc;
                g_S_LpM = wS*S_disc - (wL*L_disc + wM*M_disc);
                G_DKL = [g_LmM, g_S_LpM];
                chanNames_DKL = {'L_minus_M', 'S_minus_LpM'};
                
                [~, bestChan_DKL, DKL_b0, DKL_betas, R2_best_trial_DKL, R2_all_trial_DKL] = ...
                    compute_space_stats_trial(G_DKL, y, chanNames_DKL);
                
                for ii=1:numel(chanNames_DKL)
                    switch chanNames_DKL{ii}
                        case 'L_minus_M', DKL_bLmM=DKL_betas(ii);
                        case 'S_minus_LpM', DKL_bSLpM=DKL_betas(ii);
                    end
                end
            end
            
            fprintf('  [%s] %s | Trial LMS best=%s R2=%.3f | DKL best=%s R2=%.3f\n', ...
                sessionFolderName, unitLabel, bestChan_LMS, R2_all_trial_LMS, bestChan_DKL, R2_all_trial_DKL);
                
            allRows = [allRows; struct(...
                'date', sessionFolderName, 'unit', unitLabel, 'unit_type', unitType, 'phy_id', phyID, ...
                'LMS_best_channel', bestChan_LMS, 'LMS_R2_best_trial', R2_best_trial_LMS, 'LMS_R2_all_trial', R2_all_trial_LMS, ...
                'LMS_b0', LMS_b0, 'LMS_bL', LMS_bL, 'LMS_bM', LMS_bM, 'LMS_bS', LMS_bS, ...
                'DKL_best_channel', bestChan_DKL, 'DKL_R2_best_trial', R2_best_trial_DKL, 'DKL_R2_all_trial', R2_all_trial_DKL, ...
                'DKL_b0', DKL_b0, 'DKL_bLmM', DKL_bLmM, 'DKL_bSLpM', DKL_bSLpM)];
        end
        
        % --- SAVE PER SESSION ---
        if ~isempty(allRows)
            T_save = struct2table(allRows);
            csvName = sprintf('%s_STA_LMS_DKL_summary_trialfit.csv', sessionFolderName);
            writetable(T_save, fullfile(outDirCentral, csvName));
            writetable(T_save, fullfile(outDirLocal, csvName));
            allRowsGlobal = [allRowsGlobal; allRows];
        end
    end
    
    % --- SAVE GLOBAL SUMMARY ---
    if ~isempty(allRowsGlobal)
        T_all = struct2table(allRowsGlobal);
        globalCSV = fullfile(outRootGlobal, 'ALL_STA_LMS_DKL_summary_trialfit.csv');
        writetable(T_all, globalCSV);
        fprintf('\n✅ Wrote Global Trial Summary to %s\n', globalCSV);
    end
end

function [bestIdx, bestName, tf_intercept, tf_betas, R2_best_trial, R2_all_trial] = ...
    compute_space_stats_trial(G_trial, y, chanNames)
% compute_space_stats_trial - Fits models directly to noisy trial data.
% Operates on raw spike counts to capture trial-to-trial variance.

    nChan   = size(G_trial,2);
    bestIdx = NaN; bestName = '';
    R2_best_trial = NaN; R2_all_trial = NaN;
    tf_intercept = NaN; tf_betas = nan(1,nChan);
    
    if isempty(G_trial) || isempty(y), return; end
    
    % R2 for each channel individually
    R2_each = nan(1,nChan);
    for c = 1:nChan
        gs = double(G_trial(:,c));
        ok = isfinite(y) & isfinite(gs);
        if nnz(ok) < 10, continue; end
        
        Xs = [ones(nnz(ok),1) gs(ok)];
        bs = Xs \ y(ok);
        yhat = Xs * bs;
        
        SST = sum((y(ok) - mean(y(ok))).^2);
        SSE = sum((y(ok) - yhat).^2);
        R2_each(c) = 1 - SSE / max(SST, eps);
    end
    
    [~, idxMax] = max(R2_each);
    if ~isempty(idxMax) && isfinite(R2_each(idxMax))
        bestIdx = idxMax;
        bestName = chanNames{bestIdx};
        R2_best_trial = R2_each(bestIdx);
    end
    
    % R2 for full space (all channels)
    ok_all = isfinite(y) & all(isfinite(G_trial),2);
    if nnz(ok_all) >= 10
        Xs_all = [ones(nnz(ok_all),1) G_trial(ok_all,:)];
        bs_all = Xs_all \ y(ok_all);
        tf_intercept = bs_all(1);
        tf_betas = bs_all(2:end).';
        
        yhat_all = Xs_all * bs_all;
        SST = sum((y(ok_all) - mean(y(ok_all))).^2);
        SSE = sum((y(ok_all) - yhat_all).^2);
        R2_all_trial = 1 - SSE / max(SST, eps);
    end
end