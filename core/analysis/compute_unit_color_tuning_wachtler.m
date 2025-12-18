function unitTuning = compute_unit_color_tuning_wachtler(U, COL, config)

unitTuning = struct();

if ~isfield(U,'spk') || U.nTrials == 0
    unitTuning.fit = [];
    return;
end

hueAll = U.trials.hueID;
satAll = U.trials.satID;

nTr = numel(U.spk);
hasOutliers = isfield(U,'isOutlier') && numel(U.isOutlier) == nTr;

% time windows
sigWin = [0.05 0.15];
if isfield(config,'tuning') && isfield(config.tuning,'win')
    sigWin = config.tuning.win;
elseif isfield(U,'winEarly')
    sigWin = U.winEarly;
end

baseWin = [0.25 0.4];
if isfield(config,'tuning') && isfield(config.tuning,'baselineWin')
    baseWin = config.tuning.baselineWin;
end

durSig  = diff(sigWin);
durBase = diff(baseWin);

rateSig  = nan(nTr,1);
rateBase = nan(nTr,1);

for tt = 1:nTr
    spks = U.spk{tt};
    if isempty(spks), continue; end
    rateSig(tt)  = sum(spks >= sigWin(1)  & spks < sigWin(2))  ./ durSig;
    rateBase(tt) = sum(spks >= baseWin(1) & spks < baseWin(2)) ./ durBase;
end

rate = rateSig - rateBase;

% which saturations to use
useSats = [];
if isfield(config,'tuning') && isfield(config.tuning,'useSaturIDs') ...
        && ~isempty(config.tuning.useSaturIDs)
    useSats = config.tuning.useSaturIDs(:).';
end

satTol = 1e-3;
if isfield(config,'legacy') && isfield(config.legacy,'satTol')
    satTol = config.legacy.satTol;
end

% hue list = all hues that ever appear
hueList = unique(hueAll(isfinite(hueAll)));
hueList = hueList(:);

% saturation list: either requested sats or all sats present
if isempty(useSats)
    satList = unique(satAll(isfinite(satAll))).';
else
    satList = useSats(:).';
end

nH = numel(hueList);
nS = numel(satList);

meanRate   = nan(nH,nS);
semRate    = nan(nH,nS);
nTrialsMat = zeros(nH,nS);

validBase = isfinite(hueAll) & isfinite(satAll);
if hasOutliers
    validBase = validBase & ~U.isOutlier(:);
end
if ~isempty(useSats)
    satMask = false(size(satAll));
    for s = satList
        satMask = satMask | abs(satAll - s) < satTol;
    end
    validBase = validBase & satMask;
end

for ih = 1:nH
    h = hueList(ih);
    for is = 1:nS
        sTarget = satList(is);

        idx = isfinite(hueAll) & isfinite(satAll) & ...
              hueAll == h & abs(satAll - sTarget) < satTol;

        if hasOutliers
            idx = idx & ~U.isOutlier(:);
        end

        vals = rate(idx);
        nTrialsMat(ih,is) = numel(vals);
        if isempty(vals), continue; end

        meanRate(ih,is) = mean(vals,'omitnan');
        semRate(ih,is)  = std(vals,'omitnan') ./ sqrt(numel(vals));
    end
end

baselineHz = mean(rateBase(validBase),'omitnan');

thetaDeg = spatial_transforms.hue_to_angle(hueList, COL, config);

unitTuning.thetaID    = hueList;
unitTuning.thetaDeg   = thetaDeg;
unitTuning.saturIDs   = satList;
unitTuning.meanRate   = meanRate;
unitTuning.semRate    = semRate;
unitTuning.nTrials    = nTrialsMat;
unitTuning.baselineHz = baselineHz;

% choose which saturation(s) to use for the fit
if isempty(satList)
    unitTuning.fit = [];
    return;
end

fitMask = true(size(satList));
if ~isempty(useSats)
    fitMask = false(size(satList));
    for s = useSats
        fitMask = fitMask | abs(satList - s) < satTol;
    end
end

satIdx = find(fitMask);
if isempty(satIdx)
    satIdx = 1:numel(satList);
end

[~, kLocal] = max(satList(satIdx));
bestSatIdx = satIdx(kLocal);

theta = thetaDeg(:);
resp  = meanRate(:,bestSatIdx);

validFit = isfinite(theta) & isfinite(resp);
theta = theta(validFit);
resp  = resp(validFit);

if numel(resp) < 3 || all(resp <= 0 | isnan(resp))
    unitTuning.fit = [];
    return;
end

resp(resp < 0) = 0;

unitTuning.fit = fit_circular_tuning(theta, resp, struct());

end
