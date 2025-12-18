function DiscProbe_STA_LMSLINv6_DirectFit(DATES)
% Direct fitting version - fits LMS/DKL models directly to data
% Uses stimulus-mean responses only (no per-trial fitting)

clc;

rootDir = fileparts(mfilename('fullpath'));
addpath(genpath(rootDir));

config = load_config();
MONKEY = config.monkey;

PATHS = struct();
PATHS.baseDiscProbeLocal = config.paths.base;
PATHS.baseDiscProbesCode = config.paths.code;

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
warning('off','MATLAB:rankDeficientMatrix');

timeWin   = [0 0.2];
saturIDs  = [0.33 0.66 1];
nSat      = numel(saturIDs);
N_BOOT_REPS = 1000;

% collect everything across sessions
allRowsGlobal = [];

% ------------- keep list ----------------
keepMap = load_keep_map(PATHS.baseDiscProbeLocal);


% ------------- color metadata ----------------
lmsFile = fullfile(PATHS.baseDiscProbesCode,'lms.csv');
if ~isfile(lmsFile), error('Missing lms.csv at %s', lmsFile); end
LMS_all = readmatrix(lmsFile);

colsFile = fullfile(PATHS.baseDiscProbesCode, 'DiscProbeColsUpdated.mat');
if ~isfile(colsFile), error('Missing DiscProbeColsUpdated.mat at %s', colsFile); end
load(colsFile, 'ProbeCols', 'ProbeColIDs');
nHue = max(ProbeColIDs(:,1));
LMS_satur2 = reshape(LMS_all(1:nSat*nHue,1:3), [nSat, nHue, 3]);
nColorLMS  = size(LMS_satur2,3);

% ------------- loop over sessions ----------------
for idate = 1:numel(DATES)
    filedate = DATES{idate};
    fprintf('\n=== %s: Direct LMS & DKL fitting ===\n', filedate);

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

    trialMat = fullfile(PATHS.baseDiscProbeLocal, filedate, ...
        sprintf('%s_TrialIDs_session.mat', filedate));
    if ~isfile(trialMat)
        fprintf('  [%s] TrialIDs_session.mat not found, skipping.\n', filedate);
        continue;
    end
    T_trials = load(trialMat,'T_trials');
    T_trials = T_trials.T_trials;
    nTrials  = height(T_trials);

    hueIdx = T_trials.hueIdx;
    satVal = T_trials.satVal;
    [satIdx, okSat] = ismember(satVal, saturIDs);
    validStim = okSat & hueIdx >= 1 & hueIdx <= nHue;

    trial_LMS = nan(nTrials, nColorLMS);
    for t = 1:nTrials
        if ~validStim(t), continue; end
        trial_LMS(t,:) = squeeze(LMS_satur2(satIdx(t), hueIdx(t), :));
    end

    cols = 10:(9+nTOT);
    spkcell_all = ExptTrialsDisc(:,cols);
    spkcell_all = cellfun(@(x)x(x>=timeWin(1)&x<timeWin(2)),spkcell_all,'UniformOutput',false);
    spkcounts_all = cellfun(@numel,spkcell_all);

    outDir = fullfile(PATHS.baseDiscProbeLocal,filedate,'STA_LMSchannels');
    if ~exist(outDir,'dir'), mkdir(outDir); end
    allRows = [];

    for cc = 1:nTOT
        unitType  = iff(cc <= nSU,'SU','MU');
        unitLabel = sprintf('%s%d',unitType,unitIDs(cc));
        if ~isempty(keepMap)
            key = sprintf('%s_%s',filedate,unitLabel);
            if ~isKey(keepMap,key) || keepMap(key)==0, continue; end
        end

        phyID = unitIDs(cc);
        staFileInfo = find_sta_file(PATHS,filedate,unitType,phyID);
        if isempty(staFileInfo), continue; end

        try
            [w_chan, chanNames_sta] = grand_sta_per_color(staFileInfo.fullpath);
        catch ME
            fprintf('  [%s] %s: STA load failed: %s\n',filedate,unitLabel,ME.message);
            continue;
        end

        idxL = find(strcmpi(chanNames_sta,'L'),1);
        idxM = find(strcmpi(chanNames_sta,'M'),1);
        idxS = find(strcmpi(chanNames_sta,'S'),1);
        if isempty(idxL)||isempty(idxM)
            fprintf('  [%s] %s: missing L/M weights.\n',filedate,unitLabel);
            continue;
        end

        wL = w_chan(idxL); wM = w_chan(idxM);
        hasS = ~isempty(idxS); if hasS, wS = w_chan(idxS); end

        L_disc_all = trial_LMS(:,1);
        M_disc_all = trial_LMS(:,2);
        if hasS, S_disc_all = trial_LMS(:,3); end
        y_all = spkcounts_all(:,cc);

        if hasS
            valid = validStim & isfinite(L_disc_all)&isfinite(M_disc_all)&isfinite(S_disc_all)&isfinite(y_all);
        else
            valid = validStim & isfinite(L_disc_all)&isfinite(M_disc_all)&isfinite(y_all);
        end
        if nnz(valid)<20, continue; end

        y=y_all(valid);
        L_disc=L_disc_all(valid);
        M_disc=M_disc_all(valid);
        if hasS, S_disc=S_disc_all(valid); end
        stimID = sub2ind([nSat nHue],satIdx(valid),hueIdx(valid));
        nStim  = nSat*nHue;

        gL=wL*L_disc; gM=wM*M_disc;
        if hasS
            gS=wS*S_disc;
            G_LMS=[gL,gM,gS]; chanNames_LMS={'L','M','S'};
        else
            G_LMS=[gL,gM]; chanNames_LMS={'L','M'};
        end

        [~, bestChan_LMS, ~, ~, ~, ~, LMS_b0, LMS_betas, R2_best_stim_LMS, R2_all_stim_LMS] = ...
            compute_space_stats_direct(G_LMS, y, stimID, nSat, nHue, chanNames_LMS, N_BOOT_REPS, nStim);

        LMS_bL=NaN; LMS_bM=NaN; LMS_bS=NaN;
        for ii=1:numel(chanNames_LMS)
            switch upper(chanNames_LMS{ii})
                case 'L', LMS_bL=LMS_betas(ii);
                case 'M', LMS_bM=LMS_betas(ii);
                case 'S', LMS_bS=LMS_betas(ii);
            end
        end

        bestChan_DKL=''; R2_best_stim_DKL=NaN; R2_all_stim_DKL=NaN;
        DKL_b0=NaN; DKL_betas=[NaN NaN]; DKL_bLmM=NaN; DKL_bSLpM=NaN;
        if hasS
            g_LmM=wL*L_disc-wM*M_disc;
            g_S_LpM=wS*S_disc-(wL*L_disc+wM*M_disc);
            G_DKL=[g_LmM,g_S_LpM];
            chanNames_DKL={'L_minus_M','S_minus_LpM'};
            [~, bestChan_DKL, ~, ~, ~, ~, DKL_b0, DKL_betas, R2_best_stim_DKL, R2_all_stim_DKL] = ...
                compute_space_stats_direct(G_DKL, y, stimID, nSat, nHue, chanNames_DKL, N_BOOT_REPS, nStim);
            for ii=1:numel(chanNames_DKL)
                switch chanNames_DKL{ii}
                    case 'L_minus_M', DKL_bLmM=DKL_betas(ii);
                    case 'S_minus_LpM', DKL_bSLpM=DKL_betas(ii);
                end
            end
        else
            fprintf('  [%s] %s: no S channel → skipping DKL.\n', filedate, unitLabel);
        end

        fprintf('  [%s] %s | LMS best=%s R2_best_stim=%.3f R2_all_stim=%.3f | DKL best=%s R2_best_stim=%.3f R2_all_stim=%.3f\n', ...
            filedate, unitLabel, bestChan_LMS, R2_best_stim_LMS, R2_all_stim_LMS, bestChan_DKL, R2_best_stim_DKL, R2_all_stim_DKL);

        allRows = [allRows; struct('date',filedate,'unit',unitLabel,'unit_type',unitType,'phy_id',unitIDs(cc), ...
            'LMS_best_channel',bestChan_LMS,'LMS_R2_best_stim',R2_best_stim_LMS,'LMS_R2_all_stim',R2_all_stim_LMS, ...
            'LMS_b0',LMS_b0,'LMS_bL',LMS_bL,'LMS_bM',LMS_bM,'LMS_bS',LMS_bS, ...
            'DKL_best_channel',bestChan_DKL,'DKL_R2_best_stim',R2_best_stim_DKL,'DKL_R2_all_stim',R2_all_stim_DKL, ...
            'DKL_b0',DKL_b0,'DKL_bLmM',DKL_bLmM,'DKL_bSLpM',DKL_bSLpM)];
    end

    if ~isempty(allRows)
        T=struct2table(allRows);
        outCSV=fullfile(outDir,sprintf('%s_STA_LMS_DKL_summary_directfit.csv',filedate));
        writetable(T,outCSV);
        fprintf('  [%s] wrote %d rows to %s\n',filedate,height(T),outCSV);
        allRowsGlobal=[allRowsGlobal;allRows];
    end
end

if ~isempty(allRowsGlobal)
    T_all=struct2table(allRowsGlobal);
    allRoot=fullfile(PATHS.baseDiscProbeLocal,'STA_LMSchannels');
    if ~exist(allRoot,'dir'), mkdir(allRoot); end
    allCSV=fullfile(allRoot,'ALL_STA_LMS_DKL_summary_directfit.csv');
    writetable(T_all,allCSV);
    fprintf('\nWrote global summary with %d rows to %s\n',height(T_all),allCSV);
end

fprintf('\nDone.\n');
end

% ---------------- helpers ----------------
function val = iff(cond,a,b), if cond, val=a; else, val=b; end, end
function expFile = locate_expTrialsDisc(PATHS,MONKEY,filedate)
cands={fullfile(PATHS.baseDiscProbeLocal,filedate,sprintf('%s_%s_ExpTrialsDisc.mat',MONKEY,filedate)), ...
       fullfile(PATHS.baseDiscProbeLocal,filedate,sprintf('%s_ExpTrialsDisc.mat',filedate)), ...
       fullfile(PATHS.baseDiscProbeLocal,sprintf('%s_ExpTrialsDisc.mat',filedate))};
expFile='';
for i=1:numel(cands)
    if isfile(cands{i}), expFile=cands{i}; break; end
end
end
function info=find_sta_file(PATHS,filedate,unitType,phyID)
info=struct('fullpath','','name','');
staDir=fullfile(PATHS.baseDiscProbeLocal,filedate,'STAs');
if ~exist(staDir,'dir'),return;end
phyTag=sprintf('phy%05d',phyID);
umTag=tern(strcmpi(unitType,'SU'),'SU','MU');
d=dir(fullfile(staDir,sprintf('*%s_%s*_STA.mat',phyTag,umTag)));
if isempty(d),d=dir(fullfile(staDir,sprintf('*%s*_STA.mat',phyTag)));end
if isempty(d),return;end
info.fullpath=fullfile(staDir,d(1).name);
info.name=d(1).name;
end
function out=tern(cond,a,b),if cond,out=a;else,out=b;end,end

function [w_chan, chanNames] = grand_sta_per_color(staPath)
S = load(staPath);
if ~isfield(S,'STA'), error('STA variable not found in %s', staPath); end
STA = S.STA;
if numel(STA)>1, STA = STA(1); end
if ~isfield(STA,'sta_flat') || ~isfield(STA,'color_channels')
    error('STA.sta_flat or color_channels missing.');
end

sta_flat = STA.sta_flat;
chanNames = cellstr(STA.color_channels);
nColor = numel(chanNames);
[~, nCols] = size(sta_flat);

% safer handling to avoid reshape error
if nCols == nColor
    w_chan = mean(sta_flat,1);
elseif mod(nCols,nColor)==0
    nPix = nCols/nColor;
    sta_resh = reshape(sta_flat,[],nPix,nColor);
    w_chan = squeeze(sum(sum(sta_resh,1),2));
    w_chan = w_chan(:)';
else
    warning('STA format unexpected (%d cols, %d colors). Using column means.', nCols, nColor);
    w_chan = mean(sta_flat,1);
end

names_upper = upper(strtrim(chanNames(:)));
isRGB = (nColor==3)&&all(ismember(names_upper,{'R','G','B'}));
if isRGB
    [tL,tM,tS]=get_gun_LMS_weights();
    order=cellfun(@(n)find(strcmp({'R','G','B'},n)),names_upper);
    w_chan=[dot(w_chan,tL(order)),dot(w_chan,tM(order)),dot(w_chan,tS(order))];
    chanNames={'L','M','S'};
end
end

function [tL,tM,tS]=get_gun_LMS_weights()
rx=0.6280;ry=0.3310;gx=0.3059;gy=0.5826;bx=0.1639;by=0.0617;
Wr=18.6469/2;Wg=75.8449/2;Wb=10.5313/2;
R=cie2lms(rx,ry);G=cie2lms(gx,gy);B=cie2lms(bx,by);
M=[R;G;B]'*diag([Wr Wg Wb]);
tL=M(1,:)';tM=M(2,:)';tS=M(3,:)';
end

function v=cie2lms(a,b)
x=a;y=b;z=1-x-y;
M=[.15514 .54316 -.03286;-.15514 .45684 .03286;0 0 .01608];
lms=(M*[x;y;z])';
den=lms(1)+lms(2);
v=[lms(1)/den,lms(2)/den,lms(3)/den];
end

function [bestIdx, bestName, R2_best_trial, R2_all_trial, ...
          R2_all_boot_mean, R2_all_boot_ci, ...
          tf_intercept, tf_betas, ...
          R2_best_stim, R2_all_stim] = ...
    compute_space_stats_direct(G_trial, y, stimID, nSat, nHue, ...
                               chanNames, N_BOOT_REPS, nStim)

if nargin < 7 || isempty(N_BOOT_REPS)
    N_BOOT_REPS = 1000;
end
if nargin < 8 || isempty(nStim)
    nStim = nSat * nHue;
end

nChan   = size(G_trial,2);
bestIdx = NaN; bestName = '';
R2_best_trial = NaN; R2_all_trial = NaN;
R2_best_stim  = NaN; R2_all_stim  = NaN;
R2_all_boot_mean = NaN; R2_all_boot_ci = [NaN NaN];
tf_intercept = NaN; tf_betas = nan(1,nChan);

if isempty(G_trial) || isempty(y) || isempty(stimID)
    return;
end

% mean firing per stimulus
y_stim_full = accumarray(stimID, y, [nStim 1], @mean, NaN);
G_stim_full = nan(nStim, nChan);
for c = 1:nChan
    gc = double(G_trial(:,c));
    G_stim_full(:,c) = accumarray(stimID, gc, [nStim 1], @mean, NaN);
end

goodStim = isfinite(y_stim_full) & any(isfinite(G_stim_full), 2);
if nnz(goodStim) < 5, return; end

y_stim = y_stim_full(goodStim);
G_stim = G_stim_full(goodStim, :);

% best single channel R²
R2_each_stim = nan(1,nChan);
for c = 1:nChan
    gs = G_stim(:,c);
    good_c = isfinite(y_stim) & isfinite(gs);
    if nnz(good_c) < 5, continue; end
    ys = y_stim(good_c);
    gs_good = gs(good_c);
    Xs = [ones(nnz(good_c),1) gs_good];
    bs = Xs \ ys;
    yhat = Xs * bs;
    SST = sum((ys - mean(ys)).^2);
    SSE = sum((ys - yhat).^2);
    R2_each_stim(c) = 1 - SSE / max(SST, eps);
end

[~, idxMax] = max(R2_each_stim);
if isempty(idxMax) || ~isfinite(R2_each_stim(idxMax)), return; end

bestIdx = idxMax;
bestName = chanNames{bestIdx};
R2_best_stim = R2_each_stim(bestIdx);

% multichannel R²
good_all = isfinite(y_stim) & all(isfinite(G_stim),2);
if nnz(good_all) >= 5
    ys_all = y_stim(good_all);
    Xs_all = [ones(nnz(good_all),1) G_stim(good_all,:)];
    bs_all = Xs_all \ ys_all;
    tf_intercept = bs_all(1);
    tf_betas = bs_all(2:end).';
    yhat_all = Xs_all * bs_all;
    SST = sum((ys_all - mean(ys_all)).^2);
    SSE = sum((ys_all - yhat_all).^2);
    R2_all_stim = 1 - SSE / max(SST, eps);
end

% match legacy fieldnames for compatibility
R2_best_trial = R2_best_stim;
R2_all_trial = R2_all_stim;

% bootstrap on stimulus means
idxValid = find(good_all);
nValidStim = numel(idxValid);
if nValidStim < 10, return; end

R2_boot = nan(N_BOOT_REPS,1);
for r = 1:N_BOOT_REPS
    bootSel = idxValid(randi(nValidStim, nValidStim, 1));
    boot_y = y_stim_full(bootSel);
    boot_G = G_stim_full(bootSel,:);
    good_b = isfinite(boot_y) & all(isfinite(boot_G),2);
    if nnz(good_b) < 5, continue; end
    yb = boot_y(good_b);
    Xb = [ones(nnz(good_b),1) boot_G(good_b,:)];
    try
        bb = Xb \ yb;
        yhat_b = Xb * bb;
        SST_b = sum((yb - mean(yb)).^2);
        SSE_b = sum((yb - yhat_b).^2);
        R2_boot(r) = 1 - SSE_b / max(SST_b, eps);
    catch
        R2_boot(r) = NaN;
    end
end

if sum(isfinite(R2_boot)) >= 10
    R2_all_boot_mean = mean(R2_boot,'omitnan');
    R2_all_boot_ci = prctile(R2_boot(isfinite(R2_boot)),[2.5 97.5]);
end
end
