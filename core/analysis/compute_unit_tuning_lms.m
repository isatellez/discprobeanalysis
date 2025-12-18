function tuning3D = compute_unit_tuning_lms(U, COL, config)
% Collapse trials onto unique [L-M, S, Lum] directions using early-window rates.

tuning3D = [];

if ~isfield(U,'colIdx') || ~isfield(U,'spk') || ~isfield(U,'trials')
    warning('compute_unit_tuning_lms: missing colIdx/spk/trials, skipping.');
    return;
end

if ~isfield(COL,'axis_LM') || ~isfield(COL,'axis_S') || ~isfield(COL,'axis_Lum')
    warning('compute_unit_tuning_lms: LMS axes missing in COL, skipping.');
    return;
end

% window: prefer Wachtler config, then U.winEarly, else [0 0.2]
if isfield(config,'tuning') && isfield(config.tuning,'win') && ~isempty(config.tuning.win)
    win = config.tuning.win;
elseif isfield(U,'winEarly')
    win = U.winEarly;
else
    win = [0 0.2];
end
dur = diff(win);

% per-trial firing rate in that window from U.spk
nSpkTrials = numel(U.spk);
rate_trial = nan(nSpkTrials,1);

for tt = 1:nSpkTrials
    spks = U.spk{tt};
    if isempty(spks)
        rate_trial(tt) = 0;
    else
        use = spks >= win(1) & spks < win(2);
        rate_trial(tt) = sum(use) / dur;
    end
end

% match colIdx and rate_trial lengths
idx   = U.colIdx(:);
nIdx  = numel(idx);
nRate = numel(rate_trial);
nUse  = min(nIdx, nRate);

if nUse == 0
    return;
end

idx        = idx(1:nUse);
rate_trial = rate_trial(1:nUse);

valid = ~isnan(idx) & isfinite(rate_trial);

% respect outlier flags if present (trimmed, column, same length)
if isfield(U,'isOutlier') && numel(U.isOutlier) >= nUse
    o = U.isOutlier(1:nUse);
    o = o(:);                % force column
    valid = valid & ~o;
end

valid = valid(:);           % make sure it's a column vector

idx = idx(valid);
r   = rate_trial(valid);

if isempty(idx)
    return;
end

LM  = COL.axis_LM(idx);
S   = COL.axis_S(idx);
Lum = COL.axis_Lum(idx);

X = [LM, S, Lum];

[uniqX,~,g] = unique(X,'rows','stable');
nDir = size(uniqX,1);

meanRate = nan(nDir,1);
semRate  = nan(nDir,1);
nPerDir  = nan(nDir,1);

for k = 1:nDir
    use = (g == k);
    vals = r(use);
    vals = vals(~isnan(vals));
    if isempty(vals)
        continue;
    end
    meanRate(k) = mean(vals);
    n = numel(vals);
    semRate(k)  = std(vals) ./ max(1, sqrt(n));
    nPerDir(k)  = n;
end

tuning3D.LM   = uniqX(:,1);
tuning3D.S    = uniqX(:,2);
tuning3D.Lum  = uniqX(:,3);
tuning3D.rate = meanRate;
tuning3D.sem  = semRate;
tuning3D.n    = nPerDir;
tuning3D.win  = win;
end
