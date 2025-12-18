function DiscProbe_STA_LMSLINv5(DATES)

clc;

% make sure repo is on path (same pattern as run_wachtler_analysis)
rootDir = fileparts(mfilename('fullpath'));
addpath(genpath(rootDir));

% config + paths
config = load_config();
MONKEY = config.monkey;

PATHS = struct();
PATHS.baseDiscProbeLocal = config.paths.base;
PATHS.baseDiscProbesCode = config.paths.code;

% handle DATES: [] or 'all' means "all 6-digit date folders"
if nargin < 1 || isempty(DATES) || ...
        (ischar(DATES)   && strcmpi(strtrim(DATES),'all')) || ...
        (isstring(DATES) && strcmpi(strtrim(DATES),'all'))
    d = dir(PATHS.baseDiscProbeLocal);
    isDir  = [d.isdir];
    names  = {d(isDir).name};
    names  = names(~ismember(names, {'.','..'}));
    isDate = cellfun(@(s) numel(s)==6 && all(isstrprop(s,'digit')), names);
    DATES  = sort(names(isDate));
end

if ischar(DATES) || isstring(DATES)
    DATES = cellstr(DATES);
end

warning('off','stats:glmfit:IllConditioned');
warning('off','stats:glmfit:IterationLimit');

timeWin   = [0 0.2];    % spike counting window
saturIDs  = [0.33 0.66 1];
nSat      = numel(saturIDs);

N_BOOT_REPS = 1000;

% =========================
% keep list (VSS)
% =========================
keepMap  = containers.Map;

% primary: VSS keep list
keepFileVSS = fullfile(PATHS.baseDiscProbeLocal, 'ALL_keep_vss.csv');
% optional fallback: old global keep list
keepFileAll = fullfile(PATHS.baseDiscProbeLocal, 'ALL_keep.csv');

if isfile(keepFileVSS)
    keepFile = keepFileVSS;
elseif isfile(keepFileAll)
    keepFile = keepFileAll;
    warning('VSS keep list not found, falling back to ALL_keep.csv at %s', keepFileAll);
else
    warning('No keep list found (ALL_keep_vss.csv or ALL_keep.csv). Processing ALL units.');
    keepMap = containers.Map;  % empty; later code will skip keep filtering
    % and we just return out of this block
end

if exist('keepFile','var')
    fprintf('Loading keep list from %s...\n', keepFile);
    T_keep = readtable(keepFile);

    names      = T_keep.Properties.VariableNames;
    namesLower = lower(names);

    % session/date column
    sessField = '';
    for cand = {'session','date','datestr'}
        idx = find(strcmp(namesLower, cand{1}), 1);
        if ~isempty(idx)
            sessField = names{idx};
            break;
        end
    end
    if isempty(sessField)
        error('Keep list %s is missing a session/date column.', keepFile);
    end

    % unit label column
    unitField = '';
    for cand = {'unit','unit_label','unitname'}
        idx = find(strcmp(namesLower, cand{1}), 1);
        if ~isempty(idx)
            unitField = names{idx};
            break;
        end
    end
    if isempty(unitField)
        error('Keep list %s is missing a unit column.', keepFile);
    end

    % keep flag column (prefer VSS-specific if present)
    keepField = '';
    for cand = {'keep_vss','keep','keep_flag'}
        idx = find(strcmp(namesLower, cand{1}), 1);
        if ~isempty(idx)
            keepField = names{idx};
            break;
        end
    end
    if isempty(keepField)
        error('Keep list %s is missing a keep flag column.', keepFile);
    end

    for i = 1:height(T_keep)
        uName = string(T_keep.(unitField)(i));

        sVal = T_keep.(sessField)(i);
        if isnumeric(sVal)
            sessStr = sprintf('%06d', sVal);
        else
            sessStr = char(string(sVal));
        end

        key = sprintf('%s_%s', sessStr, uName);
        keepMap(key) = T_keep.(keepField)(i);
    end
end


% % =========================
% % keep list
% % =========================
% keepFile = fullfile(PATHS.baseDiscProbeLocal, 'ALL_keep.csv');
% keepMap  = containers.Map;
% 
% if isfile(keepFile)
%     fprintf('Loading keep list from %s...\n', keepFile);
%     T_keep = readtable(keepFile);
%     for i = 1:height(T_keep)
%         uName = string(T_keep.unit(i));
% 
%         sVal  = T_keep.session(i);
%         if isnumeric(sVal)
%             sessStr = sprintf('%06d', sVal);
%         else
%             sessStr = char(string(sVal));
%         end
% 
%         key = sprintf('%s_%s', sessStr, uName);
%         keepMap(key) = T_keep.keep(i);
%     end
% else
%     warning('ALL_keep.csv not found at %s. Processing ALL units.', keepFile);
% end


% =========================
% LMS table for disc stimuli
% =========================
lmsFile = fullfile(PATHS.baseDiscProbesCode,'lms.csv');
if ~isfile(lmsFile)
    error('Missing lms.csv at %s', lmsFile);
end
LMS_all = readmatrix(lmsFile);

% disc probe color definitions
colsFile = fullfile(PATHS.baseDiscProbesCode, 'DiscProbeColsUpdated.mat');
if ~isfile(colsFile)
    error('Missing DiscProbeColsUpdated.mat at %s', colsFile);
end
load(colsFile, 'ProbeCols', 'ProbeColIDs'); %#ok<NASGU>
nHue = max(ProbeColIDs(:,1));

if size(LMS_all,1) < nSat*nHue || size(LMS_all,2) < 3
    error('lms.csv shape is inconsistent with nSat=%d, nHue=%d', nSat, nHue);
end

% reshape to [nSat × nHue × 3] for L,M,S
LMS_satur2 = reshape(LMS_all(1:nSat*nHue,1:3), [nSat, nHue, 3]);
nColorLMS  = size(LMS_satur2,3);   % should be 3

% =========================
% main loop over sessions
% =========================
for idate = 1:numel(DATES)
    filedate = DATES{idate};
    fprintf('\n=== %s: STA → LMS & DKL vs cosine tuning ===\n', filedate);

    % load trial struct
    expFile = locate_expTrialsDisc(PATHS, MONKEY, filedate);
    if ~isfile(expFile)
        fprintf('  [%s] ExpTrialsDisc file not found, skipping.\n', filedate);
        continue;
    end
    S = load(expFile);
    ExptTrialsDisc = S.ExptTrialsDisc;
    ExptInfo       = S.ExptInfo;

    nSU     = ExptInfo.nSU;
    nMU     = ExptInfo.nMU;
    nTOT    = nSU + nMU;
    unitIDs = [ExptInfo.spk_ID_SU; ExptInfo.spk_ID_MU];

    % trial info
    trialMat = fullfile(PATHS.baseDiscProbeLocal, filedate, ...
                        sprintf('%s_TrialIDs_session.mat', filedate));
    if ~isfile(trialMat)
        fprintf('  [%s] TrialIDs_session.mat not found, skipping.\n', filedate);
        continue;
    end
    T_trials = load(trialMat, 'T_trials');
    T_trials = T_trials.T_trials;
    nTrials  = height(T_trials);

    hueIdx = T_trials.hueIdx;
    satVal = T_trials.satVal;

    % which trials have valid hue/sat
    [satIdx, okSat] = ismember(satVal, saturIDs);
    validStim = okSat & hueIdx >= 1 & hueIdx <= nHue;

    % trial LMS [nTrials × 3] (stimulus LMS, not STA weights)
    trial_LMS = nan(nTrials, nColorLMS);
    for t = 1:nTrials
        if ~validStim(t)
            continue;
        end
        si = satIdx(t);
        hi = hueIdx(t);
        trial_LMS(t,:) = squeeze(LMS_satur2(si, hi, :));
    end

    % spike counts in analysis window
    cols = 10:(9+nTOT);
    spkcell_all = ExptTrialsDisc(:, cols);
    spkcell_all = cellfun(@(x) x(x>=timeWin(1) & x<timeWin(2)), ...
                          spkcell_all, 'UniformOutput', false);
    spkcounts_all = cellfun(@numel, spkcell_all);

    outDir = fullfile(PATHS.baseDiscProbeLocal, filedate, 'STA_LMSchannels');
    if ~exist(outDir,'dir'), mkdir(outDir); end

    allRows = [];

    % -------------------------
    % per-unit loop
    % -------------------------
    for cc = 1:nTOT
        unitType  = iff(cc <= nSU, 'SU', 'MU');
        unitLabel = sprintf('%s%d', unitType, unitIDs(cc));

        % keep filter
        if ~isempty(keepMap)
            key = sprintf('%s_%s', filedate, unitLabel);
            if ~isKey(keepMap, key) || keepMap(key) == 0
                continue;
            end
        end

        % STA file
        phyID = unitIDs(cc);
        staFileInfo = find_sta_file(PATHS, filedate, unitType, phyID);
        if isempty(staFileInfo)
            continue;
        end

        % STA → LMS weights (L,M and compute S if needed)
        try
            [w_chan, chanNames_sta] = grand_sta_per_color(staFileInfo.fullpath);
        catch ME
            fprintf('  [%s] %s: STA load/grand failed: %s\n', ...
                    filedate, unitLabel, ME.message);
            continue;
        end

        idxL = find(strcmpi(chanNames_sta,'L'),1); %grand weight L
        idxM = find(strcmpi(chanNames_sta,'M'),1);
        idxS = find(strcmpi(chanNames_sta,'S'),1);

        if isempty(idxL) || isempty(idxM)
            fprintf('  [%s] %s: STA did not provide L/M weights, skipping.\n', ...
                    filedate, unitLabel);
            continue;
        end
        wL = w_chan(idxL);
        wM = w_chan(idxM);
        hasS = ~isempty(idxS);
        if hasS
            wS = w_chan(idxS);
        end

        if nColorLMS < 2
            fprintf('  [%s] %s: trial LMS has <2 channels, skipping.\n', ...
                    filedate, unitLabel);
            continue;
        end

        % valid trials: good stimulus, finite LMS + spikes
        L_disc_all = trial_LMS(:,1);
        M_disc_all = trial_LMS(:,2);
        if hasS
            S_disc_all = trial_LMS(:,3);
        end
        y_all = spkcounts_all(:, cc);

        if hasS
            valid = validStim & isfinite(L_disc_all) & isfinite(M_disc_all) & ...
                    isfinite(S_disc_all) & isfinite(y_all);
        else
            valid = validStim & isfinite(L_disc_all) & isfinite(M_disc_all) & ...
                    isfinite(y_all);
        end

        if nnz(valid) < 20
            continue;
        end

        y      = y_all(valid);
        L_disc = L_disc_all(valid);
        M_disc = M_disc_all(valid);
        if hasS
            S_disc = S_disc_all(valid);
        end

        % stimulus identity per trial
        stimID = sub2ind([nSat nHue], satIdx(valid), hueIdx(valid));
        nStim  = nSat * nHue;

        % =========================
        % LMS model: channels = L, M, (S)
        % =========================
        gL = wL * L_disc;
        gM = wM * M_disc;
        if hasS
            gS = wS * S_disc;
            G_LMS_trial   = [gL, gM, gS]; %this is the generator matrix
            chanNames_LMS = {'L','M','S'};
        else
            G_LMS_trial   = [gL, gM];
            chanNames_LMS = {'L','M'};
        end

        [bestIdx_LMS, bestChan_LMS, ...
         R2_best_LMS, R2_all_LMS, ...
         R2_all_LMS_boot_mean, R2_all_LMS_boot_ci, ...
         LMS_b0, LMS_betas] = ...
            compute_space_stats(G_LMS_trial, y, stimID, nSat, nHue, ...
                                chanNames_LMS, N_BOOT_REPS, nStim);

        % map LMS betas by name
        LMS_bL = NaN; LMS_bM = NaN; LMS_bS = NaN;
        for ii = 1:numel(chanNames_LMS)
            switch upper(chanNames_LMS{ii})
                case 'L', LMS_bL = LMS_betas(ii);
                case 'M', LMS_bM = LMS_betas(ii);
                case 'S', LMS_bS = LMS_betas(ii);
            end
        end

        % =========================
        % DKL-like model: channels = L-M, S-(L+M)
        % =========================
        bestIdx_DKL          = NaN;
        bestChan_DKL         = '';
        R2_best_DKL          = NaN;
        R2_all_DKL           = NaN;
        R2_all_DKL_boot_mean = NaN;
        R2_all_DKL_boot_ci   = [NaN NaN];
        DKL_b0               = NaN;
        DKL_betas            = [NaN NaN];

        if hasS
            g_LmM   = wL * L_disc - wM * M_disc;
            g_S_LpM = wS * S_disc - (wL * L_disc + wM * M_disc);
            G_DKL_trial   = [g_LmM, g_S_LpM]; %the other generator matrix. this is already the generator signal!
            chanNames_DKL = {'L_minus_M','S_minus_LpM'};

            [bestIdx_DKL, bestChan_DKL, ...
             R2_best_DKL, R2_all_DKL, ...
             R2_all_DKL_boot_mean, R2_all_DKL_boot_ci, ...
             DKL_b0, DKL_betas] = ...
                compute_space_stats(G_DKL_trial, y, stimID, nSat, nHue, ...
                                    chanNames_DKL, N_BOOT_REPS, nStim);
        else
            fprintf('  [%s] %s: no S channel → skipping DKL analysis.\n', ...
                    filedate, unitLabel);
        end

        % map DKL betas
        DKL_bLmM  = NaN;
        DKL_bSLpM = NaN;
        if exist('chanNames_DKL','var')
            for ii = 1:numel(chanNames_DKL)
                switch chanNames_DKL{ii}
                    case 'L_minus_M',   DKL_bLmM  = DKL_betas(ii);
                    case 'S_minus_LpM', DKL_bSLpM = DKL_betas(ii);
                end
            end
        end

        fprintf(['  [%s] %s | LMS: best=%s R2_best=%.3f R2_all=%.3f ', ...
                 '| DKL: best=%s R2_best=%.3f R2_all=%.3f\n'], ...
                filedate, unitLabel, ...
                bestChan_LMS, R2_best_LMS, R2_all_LMS, ...
                bestChan_DKL, R2_best_DKL, R2_all_DKL);

        ciL = R2_all_LMS_boot_ci;
        if numel(ciL) < 2
            ciL = [NaN NaN];
        end

        ciD = R2_all_DKL_boot_ci;
        if numel(ciD) < 2
            ciD = [NaN NaN];
        end

        allRows = [allRows; struct( ...
            'date', filedate, ...
            'unit', unitLabel, ...
            'unit_type', unitType, ...
            'phy_id', unitIDs(cc), ...
            'LMS_best_channel', bestChan_LMS, ...
            'LMS_best_idx', bestIdx_LMS, ...
            'LMS_R2_best', R2_best_LMS, ...
            'LMS_R2_all', R2_all_LMS, ...
            'LMS_R2_all_boot_mean', R2_all_LMS_boot_mean, ...
            'LMS_R2_all_boot_lo', ciL(1), ...
            'LMS_R2_all_boot_hi', ciL(2), ...
            'LMS_b0', LMS_b0, ...
            'LMS_bL', LMS_bL, ...
            'LMS_bM', LMS_bM, ...
            'LMS_bS', LMS_bS, ...
            'DKL_best_channel', bestChan_DKL, ...
            'DKL_best_idx', bestIdx_DKL, ...
            'DKL_R2_best', R2_best_DKL, ...
            'DKL_R2_all', R2_all_DKL, ...
            'DKL_R2_all_boot_mean', R2_all_DKL_boot_mean, ...
            'DKL_R2_all_boot_lo', ciD(1), ...
            'DKL_R2_all_boot_hi', ciD(2), ...
            'DKL_b0', DKL_b0, ...
            'DKL_bLmM', DKL_bLmM, ...
            'DKL_bSLpM', DKL_bSLpM)];
    end  % end for cc

    if ~isempty(allRows)
        T = struct2table(allRows);
        outCSV = fullfile(outDir, sprintf('%s_STA_LMS_DKL_summary.csv', filedate));
        writetable(T, outCSV);
        fprintf('  [%s] wrote %d rows to %s\n', filedate, height(T), outCSV);
    end
end  % end for idate

fprintf('\nDone.\n');
end

% ==========================================
% helpers
% ==========================================

function val = iff(cond, a, b)
if cond, val = a; else, val = b; end
end

function expFile = locate_expTrialsDisc(PATHS, MONKEY, filedate)
cands = {
    fullfile(PATHS.baseDiscProbeLocal, filedate, ...
             sprintf('%s_%s_ExpTrialsDisc.mat', MONKEY, filedate)), ...
    fullfile(PATHS.baseDiscProbeLocal, filedate, ...
             sprintf('%s_ExpTrialsDisc.mat', filedate)), ...
    fullfile(PATHS.baseDiscProbeLocal, ...
             sprintf('%s_ExpTrialsDisc.mat', filedate))};
expFile = '';
for i = 1:numel(cands)
    if isfile(cands{i})
        expFile = cands{i};
        break;
    end
end
end

function info = find_sta_file(PATHS, filedate, unitType, phyID)
info = struct('fullpath','','name','');
staDir = fullfile(PATHS.baseDiscProbeLocal, filedate, 'STAs');
if ~exist(staDir,'dir')
    return;
end

phyTag = sprintf('phy%05d', phyID);
umTag  = tern(strcmpi(unitType,'SU'),'SU','MU');

d = dir(fullfile(staDir, sprintf('*%s_%s*_STA.mat', phyTag, umTag)));
if isempty(d)
    d = dir(fullfile(staDir, sprintf('*%s*_STA.mat', phyTag)));
end
if isempty(d)
    return;
end

info.fullpath = fullfile(staDir, d(1).name);
info.name     = d(1).name;
end

function out = tern(cond,a,b)
if cond, out = a; else, out = b; end
end

function [w_chan, chanNames] = grand_sta_per_color(staPath)
S = load(staPath);
if ~isfield(S,'STA')
    error('STA variable not found in %s', staPath);
end
STA = S.STA;
if numel(STA) > 1
    STA = STA(1);
end
if ~isfield(STA,'sta_flat') || ~isfield(STA,'color_channels')
    error('STA.sta_flat or STA.color_channels missing.');
end

sta_flat = STA.sta_flat;
[~, nCols] = size(sta_flat);

baseChanNames = cellstr(STA.color_channels);
nColor        = numel(baseChanNames);

if mod(nCols, nColor) ~= 0
    error('sta_flat columns (%d) not divisible by nColor (%d)', nCols, nColor);
end

nPix = nCols / nColor;
sta_resh = reshape(sta_flat, [], nPix, nColor);

w_rgb = squeeze(sum(sum(sta_resh, 1), 2));
w_rgb = w_rgb(:);

w_chan    = w_rgb(:)';
chanNames = baseChanNames;

names_upper = upper(strtrim(baseChanNames(:)));

has_tL = isfield(STA,'transform_L') && numel(STA.transform_L) == nColor;
has_tM = isfield(STA,'transform_M') && numel(STA.transform_M) == nColor;
has_tS = isfield(STA,'transform_S') && numel(STA.transform_S) == nColor;

tL = []; tM = []; tS = [];

if has_tL && has_tM
    tL = STA.transform_L(:);
    tM = STA.transform_M(:);
end

isRGB = (nColor == 3) && all(ismember(names_upper, {'R','G','B'}));
if isRGB
    [tL_cal, tM_cal, tS_cal] = get_gun_LMS_weights();
    order = zeros(3,1);
    for k = 1:3
        switch names_upper{k}
            case 'R', order(k) = 1;
            case 'G', order(k) = 2;
            case 'B', order(k) = 3;
        end
    end
    tL_cal = tL_cal(order);
    tM_cal = tM_cal(order);
    tS_cal = tS_cal(order);

    if isempty(tL)
        tL = tL_cal;
    end
    if isempty(tM)
        tM = tM_cal;
    end
    if ~has_tS
        tS = tS_cal;
    else
        tS = STA.transform_S(:);
    end
elseif has_tS
    tS = STA.transform_S(:);
end

if ~isempty(tL) && ~isempty(tM)
    wL = dot(w_rgb, tL);
    wM = dot(w_rgb, tM);
    if ~isempty(tS)
        wS = dot(w_rgb, tS);
        w_chan    = [wL, wM, wS];
        chanNames = {'L','M','S'};
    else
        w_chan    = [wL, wM];
        chanNames = {'L','M'};
    end
end
end

function [tL, tM, tS] = get_gun_LMS_weights()
rx=0.6280; ry=0.3310; gx=0.3059; gy=0.5826; bx=0.1639; by=0.0617;
Wr=18.6469/2; Wg=75.8449/2; Wb=10.5313/2;

R = cie2lms(rx,ry);
G = cie2lms(gx,gy);
B = cie2lms(bx,by);

Lr=R(1); Lg=G(1); Lb=B(1);
Mr=R(2); Mg=G(2); Mb=B(2);
Sr=R(3); Sg=G(3); Sb=B(3);

M_rgb2lms = [Lr Lg Lb; Mr Mg Mb; Sr Sg Sb] * diag([Wr Wg Wb]);

tL = M_rgb2lms(1,:).';
tM = M_rgb2lms(2,:).';
tS = M_rgb2lms(3,:).';
end

function v = cie2lms(a,b)
x=a; y=b; z=1-x-y;
M=[ .15514 .54316 -.03286; -.15514 .45684 .03286; 0 0 .01608];
lms = (M*[x;y;z]).';
den = lms(1)+lms(2);
v = [lms(1)/den, lms(2)/den, lms(3)/den];
end

% ==========================================
% core analysis helper for one "space"
% ==========================================
function [bestIdx, bestName, R2_best, R2_all, ...
          R2_all_boot_mean, R2_all_boot_ci, ...
          tf_intercept, tf_betas] = ...
    compute_space_stats(G_trial, y, stimID, nSat, nHue, ...
                        chanNames, N_BOOT_REPS, nStim)

if nargin < 8 || isempty(N_BOOT_REPS)
    N_BOOT_REPS = 1000;
end
if nargin < 9
    nStim = nSat * nHue;
end

nChan   = size(G_trial,2);
bestIdx = NaN;
bestName = '';
R2_best = NaN;
R2_all  = NaN;
R2_all_boot_mean = NaN;
R2_all_boot_ci   = [NaN NaN];

tf_intercept = NaN;
tf_betas     = nan(1, nChan);

if isempty(G_trial) || isempty(y)
    return;
end

% tuning-level means from data
y_tune_raw = accumarray(stimID, y, [nStim 1], @mean, NaN);

% cosine fit on tuning
y_tune = fit_cosine_from_tuning(y_tune_raw, nSat, nHue);
if isempty(y_tune)
    return;
end

% mean generators per stimulus
G_tune = nan(nStim, nChan);
for c = 1:nChan
    vals_c = double(G_trial(:,c));
    G_tune(:,c) = accumarray(stimID, vals_c, [nStim 1], @mean, NaN);
end

% best single-channel R² (still computed, using cosine-tuning y_tune)
R2_each = nan(1, nChan);
for c = 1:nChan
    g_c = G_tune(:,c);
    good = isfinite(y_tune) & isfinite(g_c);
    if nnz(good) < 5
        continue;
    end
    yc = y_tune(good);
    gc = g_c(good);
    Xc = [ones(nnz(good),1) gc];
    b  = Xc \ yc;
    yhat = Xc * b;
    SST = sum((yc - mean(yc)).^2);
    SSE = sum((yc - yhat).^2);
    R2_each(c) = 1 - SSE / max(SST, eps);
end

[~, idxMax] = max(R2_each);
if ~isempty(idxMax) && isfinite(R2_each(idxMax))
    bestIdx  = idxMax;
    R2_best  = R2_each(idxMax);
    bestName = chanNames{bestIdx};
else
    return;
end

% pooled multi-channel R² and multivariate betas (this is your TF)
good_full = isfinite(y_tune) & all(isfinite(G_tune),2);
if nnz(good_full) >= 5
    y_full = y_tune(good_full);
    X_full = [ones(nnz(good_full),1) G_tune(good_full,:)];
    b_full = X_full \ y_full;     % length 1 + nChan
    tf_intercept = b_full(1);
    tf_betas     = b_full(2:end).';

    yhat_full = X_full * b_full;
    SST_full  = sum((y_full - mean(y_full)).^2);
    SSE_full  = sum((y_full - yhat_full).^2);
    R2_all    = 1 - SSE_full / max(SST_full, eps);
end

% bootstrap CIs on pooled R²
nValidTrials = numel(y);
if nValidTrials < 10
    return;
end

R2_boot = nan(N_BOOT_REPS,1);
for r = 1:N_BOOT_REPS
    bootIdx = randi(nValidTrials, nValidTrials, 1);

    bootStimID  = stimID(bootIdx);
    boot_y      = y(bootIdx);
    boot_Gtrial = G_trial(bootIdx,:);

    y_tune_raw_b = accumarray(bootStimID, boot_y, [nStim 1], @mean, NaN);
    y_tune_b     = fit_cosine_from_tuning(y_tune_raw_b, nSat, nHue);
    if isempty(y_tune_b)
        continue;
    end

    G_tune_b = nan(nStim, nChan);
    for c = 1:nChan
        vals_b = double(boot_Gtrial(:,c));
        G_tune_b(:,c) = accumarray(bootStimID, vals_b, [nStim 1], @mean, NaN);
    end

    good_b = isfinite(y_tune_b) & all(isfinite(G_tune_b),2);
    if nnz(good_b) < 5
        continue;
    end
    yb = y_tune_b(good_b);
    Xb = [ones(nnz(good_b),1) G_tune_b(good_b,:)];
    try
        b_b   = Xb \ yb;
        yhatb = Xb * b_b;
        SST_b = sum((yb - mean(yb)).^2);
        SSE_b = sum((yb - yhatb).^2);
        R2_boot(r) = 1 - SSE_b / max(SST_b, eps);
    catch
        R2_boot(r) = NaN;
    end
end

if sum(isfinite(R2_boot)) >= 10
    R2_all_boot_mean = mean(R2_boot, 'omitnan');
    R2_all_boot_ci   = prctile(R2_boot(isfinite(R2_boot)), [2.5 97.5]);
end
end

% ==========================================
% cosine fit helper (uses tuning means)
% ==========================================
function y_cos = fit_cosine_from_tuning(y_tune_raw, nSat, nHue)
nStim = numel(y_tune_raw);
if nStim ~= nSat * nHue
    y_cos = [];
    return;
end

y_cos = nan(size(y_tune_raw));

[~, HueGrid] = ind2sub([nSat nHue], (1:nStim).');
theta = 2*pi*(HueGrid-1)/nHue;

good = isfinite(y_tune_raw);
if nnz(good) < 5
    y_cos = [];
    return;
end

th_g = theta(good);
y_g  = y_tune_raw(good);

X = [ones(nnz(good),1), cos(th_g), sin(th_g)];
b = X \ y_g;

X_all = [ones(nStim,1), cos(theta), sin(theta)];
y_cos = X_all * b;
end
