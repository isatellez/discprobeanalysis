function [peakClassStr, isTunedFourier, kMax, fracMax] = classify_vss_fourier_peak(P, FT, config)
peakClassStr   = "flat";
isTunedFourier = false;
kMax           = NaN;
fracMax        = NaN;

if isempty(FT) || ~isfield(FT,'over') || isempty(FT.over)
    return;
end

F = FT.over;

if ~isfield(F,'varExp') || isempty(F.varExp)
    return;
end

K = numel(F.varExp);
kMax = F.k_max_var;
if isnan(kMax) || kMax < 1 || kMax > K
    return;
end

fracMax = F.varExp(kMax);

% thresholds
alpha    = 0.05;
minFrac  = 0.20;

if isfield(config,'vss')
    if isfield(config.vss,'alpha'),   alpha   = config.vss.alpha;   end
    if isfield(config.vss,'minFrac'), minFrac = config.vss.minFrac; end
end

% significance from peak_model stats (now using parametric P.p_tuned)
p_tuned   = NaN;
p_bimodal = NaN;
if isfield(P,'p_tuned'),   p_tuned   = P.p_tuned;   end
if isfield(P,'p_bimodal'), p_bimodal = P.p_bimodal; end

sigTuned = ~isnan(p_tuned)   && (p_tuned   < alpha);
sigBi    = ~isnan(p_bimodal) && (p_bimodal < alpha);

if ~sigTuned
    return;
end

if kMax > 2 || fracMax < minFrac
    return;
end

if kMax == 1
    peakClassStr   = "unimodal";
    isTunedFourier = true;
elseif kMax == 2
    if sigBi || isnan(p_bimodal)
        peakClassStr   = "bimodal";
        isTunedFourier = true;
    else
        peakClassStr   = "flat";
        isTunedFourier = false;
    end
end

end
