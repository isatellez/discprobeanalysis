function stats = compute_hue_means(U)
% mean firing per hue (early window, across all sats)

hue = U.trials.hueID;
hue(isnan(hue)) = [];

hues = unique(hue);
hues = hues(:);

nH = numel(hues);
rate_mean = nan(nH,1);
rate_sem  = nan(nH,1);
n_trials  = nan(nH,1);

win = U.winEarly;
dur = diff(win);

hasOutliers = isfield(U, 'isOutlier') && numel(U.isOutlier) == height(U.trials);

for ih = 1:nH
    h = hues(ih);
    use = U.trials.hueID == h;
    if hasOutliers
        use = use & ~U.isOutlier(:);
    end

    idx = find(use);
    if isempty(idx)
        continue;
    end

    r = nan(numel(idx),1);
    for ii = 1:numel(idx)
        t = idx(ii);
        spks = U.spk{t};
        if isempty(spks)
            r(ii) = 0;
        else
            r(ii) = sum(spks >= win(1) & spks < win(2)) ./ dur;
        end
    end

    rate_mean(ih) = mean(r, 'omitnan');
    rate_sem(ih)  = std(r, 'omitnan') ./ sqrt(numel(r));
    n_trials(ih)  = numel(r);
end

stats.hues      = hues;
stats.rate_mean = rate_mean;
stats.rate_sem  = rate_sem;
stats.n_trials  = n_trials;
stats.win       = win;

end
