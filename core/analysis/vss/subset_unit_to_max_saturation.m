function U_high = subset_unit_to_max_saturation(U, maxSatVal)
U_high = U;

if nargin < 2 || isempty(maxSatVal)
    satAll = U.trials.satID;
    satAll = satAll(~isnan(satAll));
    if isempty(satAll)
        U_high.nTrials = 0;
        return;
    end
    maxSatVal = max(satAll);
end

satID = U.trials.satID;
satTol = 1e-3;
use = abs(satID - maxSatVal) < satTol;

if isfield(U, 'isOutlier') && numel(U.isOutlier) == height(U.trials)
    use = use & ~U.isOutlier(:);
end

idx = find(use);
if isempty(idx)
    U_high.trials  = U.trials([]);
    U_high.spk     = U.spk([]);
    U_high.nTrials = 0;
    if isfield(U_high, 'isOutlier')
        U_high.isOutlier = U_high.isOutlier([]);
    end
    return;
end

U_high.trials  = U.trials(idx, :);
U_high.spk     = U.spk(idx);
U_high.nTrials = numel(idx);

if isfield(U, 'isOutlier')
    U_high.isOutlier = U.isOutlier(idx);
end
end
