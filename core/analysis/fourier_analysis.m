function FT = fourier_analysis(U, config)
% Fourier analysis of hue tuning (overall and per-saturation).

FT = struct();
if U.nTrials == 0
    return;
end

% defaults from original script
FOURIER.maxHarmonics = 8;
FOURIER.useHann      = true;
FOURIER.detrend      = true;

if isfield(config, 'fourier')
    fcfg = config.fourier;
    if isfield(fcfg,'maxHarmonics'), FOURIER.maxHarmonics = fcfg.maxHarmonics; end
    if isfield(fcfg,'useHann'),      FOURIER.useHann      = logical(fcfg.useHann); end
    if isfield(fcfg,'detrend'),      FOURIER.detrend      = logical(fcfg.detrend); end
end

hueID = U.trials.hueID;
satID = U.trials.satID;

% use early window to match other tuning stats
win = U.winEarly;
dur = diff(win);

nTr = numel(U.spk);
rateEarly = nan(nTr,1);
for tt = 1:nTr
    spks = U.spk{tt};
    if isempty(spks)
        rateEarly(tt) = 0;
    else
        spks = spks(spks >= win(1) & spks < win(2));
        rateEarly(tt) = numel(spks) / dur;
    end
end

% overall mean per hue
hueVals = unique(hueID(~isnan(hueID)));
hueVals = sort(hueVals(:).');
nHue    = numel(hueVals);

R_h = nan(nHue,1);
for i = 1:nHue
    h = hueVals(i);
    idx = (hueID == h);
    R_h(i) = mean(rateEarly(idx), 'omitnan');
end

% per-sat, per-hue matrix R_hs(:,si)
satVals = unique(satID(~isnan(satID)));
satVals = sort(satVals(:).');
nS      = numel(satVals);

R_hs = nan(nHue, nS);
for si = 1:nS
    s = satVals(si);
    for hi = 1:nHue
        h = hueVals(hi);
        idx = (hueID == h) & (satID == s);
        R_hs(hi,si) = mean(rateEarly(idx), 'omitnan');
    end
end

% overall Fourier over hues
FT_over = fourier_over_hues_local(R_h, FOURIER.maxHarmonics, FOURIER.useHann, FOURIER.detrend);

% per-sat Fourier
FT_cells = cell(nS,1);
for si = 1:nS
    yi = R_hs(:,si);
    if any(~isfinite(yi)), yi(~isfinite(yi)) = 0; end
    FT_cells{si} = fourier_over_hues_local(yi, FOURIER.maxHarmonics, FOURIER.useHann, FOURIER.detrend);
end

% make proper struct array
names = fieldnames(FT_cells{1});
FT_sat = repmat(struct(), nS, 1);
for si = 1:nS
    for f = 1:numel(names)
        nm = names{f};
        FT_sat(si).(nm) = FT_cells{si}.(nm);
    end
end

% derived overall metrics (mirroring your script)
K = min(FOURIER.maxHarmonics, floor(nHue/2));
FT_over.amp_norm  = FT_over.amp ./ max(FT_over.DC, eps);
FT_over.frac_noDC = FT_over.varExp;

w = FOURIER.useHann * hann(nHue,'periodic') + (~FOURIER.useHann) * ones(nHue,1);
PX   = abs(fft(R_h(:).*w)).^2;
totP = sum(PX);
noDC = sum(PX(2:end));

FT_over.DC_powerFrac      = PX(1)/totP;
FT_over.allPowerFrac_noDC = noDC/totP;
FT_over.frac_total        = (PX(2:1+K)/totP).';
FT_over.frac12_noDC       = sum(FT_over.frac_noDC(1:min(2,end)));
FT_over.frac12_wDC        = sum(PX(2:1+min(2,K)))/totP;
FT_over.fracTotal_wDC     = (PX(1) + sum(PX(2:1+K)))/totP;
[~, FT_over.k_max_amp]    = max(FT_over.amp);
[~, FT_over.k_max_var]    = max(FT_over.varExp);
FT_over.pref_deg_h1       = mod(-FT_over.phase_deg(1), 360);
FT_over.axis_deg_h2       = NaN;
if numel(FT_over.phase_deg) >= 2
    FT_over.axis_deg_h2 = mod(-FT_over.phase_deg(2)/2, 180);
end

FT.over      = FT_over;
FT.perSat    = FT_sat;
FT.saturIDs  = satVals;
FT.hueVals   = hueVals;

end

function FT = fourier_over_hues_local(y, maxHarm, useHann, doDetrend)

y = y(:);
n = numel(y);

if doDetrend
    y0 = y - mean(y, 'omitnan');
else
    y0 = y;
end
y0(~isfinite(y0)) = 0;

if useHann
    w = hann(n,'periodic');
else
    w = ones(n,1);
end

X = fft(y0 .* w);
DC = real(X(1))/n;

K  = min(maxHarm, floor(n/2));
harmonics = (1:K).';

ck = X(2:1+K);           % skip DC
amp = 2*abs(ck)/n;
ph  = angle(ck);
phase_deg = rad2deg(ph);

power_k = abs(ck).^2;
totPower_noDC = sum(power_k);
if totPower_noDC > 0
    varExp = power_k / totPower_noDC;
else
    varExp = zeros(size(power_k));
end

FT = struct();
FT.harmonics = harmonics;
FT.amp       = amp(:).';
FT.phase_deg = phase_deg(:).';
FT.varExp    = varExp(:).';
FT.DC        = DC;

end
