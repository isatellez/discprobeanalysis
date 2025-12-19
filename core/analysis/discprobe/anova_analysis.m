function A = anova_analysis(U, config)
% two-way ANOVA: hue × saturation on early window rate

win = U.winEarly;
dur = diff(win);

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

hueID = U.trials.hueID;
satID = U.trials.satID;

hasOutliers = isfield(U, 'isOutlier') && numel(U.isOutlier) == nTr;

valid = isfinite(rate) & isfinite(hueID) & isfinite(satID);
if hasOutliers
    valid = valid & ~U.isOutlier(:);
end

rate  = rate(valid);
hueID = hueID(valid);
satID = satID(valid);

if numel(rate) < 10
    A = struct();
    A.p_hue = NaN;
    A.p_sat = NaN;
    A.p_int = NaN;
    A.tbl   = {};
    A.n     = numel(rate);
    return;
end

[p, tbl, statsAn] = anovan(rate, {hueID, satID}, ...
    'model', 'full', ...
    'varnames', {'hue','sat'}, ...
    'display', 'off');

A.p_hue  = p(1);
A.p_sat  = p(2);
A.p_int  = p(3);
A.tbl    = tbl;
A.stats  = statsAn;
A.n      = numel(rate);

end
