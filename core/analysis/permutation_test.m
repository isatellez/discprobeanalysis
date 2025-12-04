function PT = permutation_test(U, config)
% permutation test on hue tuning structure
% null: shuffle hue labels across trials

if ~isfield(U, 'winEarly')
    error('U.winEarly missing.');
end

win = U.winEarly;
dur = diff(win);

nPerm = 1000;
alpha = 0.05;
if isfield(config, 'pt')
    if isfield(config.pt, 'nPerm'), nPerm = config.pt.nPerm; end
    if isfield(config.pt, 'alpha'), alpha = config.pt.alpha; end
end

hueID = U.trials.hueID;
hasOutliers = isfield(U, 'isOutlier') && numel(U.isOutlier) == height(U.trials);

nTr = numel(U.spk);
rate = nan(nTr,1);

for tt = 1:nTr
    spks = U.spk{tt};
    if isempty(spks)
        rate(tt) = 0;
    else
        rate(tt) = sum(spks >= win(1) & spks < win(2)) ./ dur;
    end
end

valid = isfinite(rate) & ~isnan(hueID);
if hasOutliers
    valid = valid & ~U.isOutlier(:);
end

rate = rate(valid);
hueID = hueID(valid);

hues = unique(hueID);
hues = hues(:);
nH   = numel(hues);

% observed tuning
obsMean = nan(nH,1);
for ih = 1:nH
    h = hues(ih);
    obsMean(ih) = mean(rate(hueID == h), 'omitnan');
end

obsStat = var(obsMean, 'omitnan');

% null distribution
nullStat = nan(nPerm,1);
nV = numel(rate);

for p = 1:nPerm
    permHue = hueID(randperm(nV));
    mPerm = nan(nH,1);
    for ih = 1:nH
        h = hues(ih);
        mPerm(ih) = mean(rate(permHue == h), 'omitnan');
    end
    nullStat(p) = var(mPerm, 'omitnan');
end

pVal = mean(nullStat >= obsStat);

PT.hues      = hues;
PT.obsMean   = obsMean;
PT.obsStat   = obsStat;
PT.nullStat  = nullStat;
PT.p         = pVal;
PT.alpha     = alpha;
PT.nPerm     = nPerm;
PT.winEarly  = win;

end
