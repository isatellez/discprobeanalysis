function U = detect_outliers(U, config)
% mark and drop trial outliers based on early window firing rates.

if ~isfield(U, 'spk') || ~isfield(U, 'trials')
    return;
end

hueID = U.trials.hueID;
satID = U.trials.satID;

win = U.winEarly;
dur = diff(win); %duration 

nTr = numel(U.spk);
rateEarly = nan(nTr,1);

for tt = 1:nTr
    spks = U.spk{tt};
    if isempty(spks)
        rateEarly(tt) = 0;
    else
        spks = spks(spks >= win(1) & spks < win(2));
        rateEarly(tt) = numel(spks) / dur; %average spikes
    end
end

% default is tukey 
useNorm = false;
method  = 'tukey';   % 'tukey' or 'mad'
minN    = 5; %min number of spikes in window to do analysis 

if isfield(config, 'boxPlots')
    if isfield(config.boxPlots, 'useNorm')
        useNorm = logical(config.boxPlots.useNorm);
    end
    if isfield(config.boxPlots, 'method')
        method = lower(config.boxPlots.method);
    end
    if isfield(config.boxPlots, 'minN')
        minN = config.boxPlots.minN;
    end
end

% normalization of firing rates if normalization was chosen 
if useNorm
    r = rateEarly;
    r = r - min(r);
    denom = max(r);
    if denom <= 0
        vals = zeros(size(r));
    else
        vals = r ./ denom;
    end
else
    vals = rateEarly;
end

is_out = false(nTr,1); %assume everyone is in at the beginning

satVals = unique(satID(~isnan(satID))); %unique sats
hueVals = unique(hueID(~isnan(hueID))); %unique hue values

for si = 1:numel(satVals)
    s = satVals(si);
    for hi = 1:numel(hueVals)
        h = hueVals(hi);

        idx = find(satID == s & hueID == h & isfinite(vals));
        if numel(idx) < minN
            continue;
        end

        x = vals(idx);

        switch method
            case 'tukey'
                q = quantile(x, [0.25 0.75]);
                I = q(2) - q(1);
                if I <= 0
                    bad = false(size(idx));
                else
                    lo = q(1) - 1.5*I;
                    hi = q(2) + 1.5*I;
                    bad = (x < lo) | (x > hi);
                end

            case 'mad'
                med = median(x);
                if exist('mad', 'file')
                    srob = 1.4826 * mad(x, 1);
                else
                    srob = 1.4826 * median(abs(x - med));
                end
                if srob <= 0
                    bad = false(size(idx));
                else
                    z = abs(x - med) / srob;
                    bad = z > 3.5;
                end

            otherwise
                bad = false(size(idx));
        end

        if any(bad)
            is_out(idx(bad)) = true;
        end
    end
end

keep = ~is_out;

if ~any(is_out)
    U.isOutlier = is_out;
    return;
end

U.isOutlier = is_out;

U.spk     = U.spk(keep);
if isfield(U, 'mrk') && numel(U.mrk) == nTr
    U.mrk = U.mrk(keep);
end

U.trials  = U.trials(keep,:);
U.nTrials = sum(keep);

end
