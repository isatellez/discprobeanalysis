function FT = fourier_analysis(U, config)
% Discrete Fourier analysis of hue tuning (early window)

hm = compute_hue_means(U);
hues = hm.hues;
resp = hm.rate_mean;

if all(~isfinite(resp)) || all(resp == 0 | isnan(resp))
    FT = struct();
    FT.hues       = hues;
    FT.resp       = resp;
    FT.harmonics  = [];
    FT.amp        = [];
    FT.phase_deg  = [];
    FT.powerFrac  = [];
    return;
end

resp = resp(:);
nH   = numel(resp);

useHann  = true;
detrendY = true;
maxHarm  = min(8, floor(nH/2));

if isfield(config, 'fourier')
    if isfield(config.fourier, 'useHann'),  useHann  = config.fourier.useHann; end
    if isfield(config.fourier, 'detrend'),  detrendY = config.fourier.detrend; end
    if isfield(config.fourier, 'maxHarmonics')
        maxHarm = min(config.fourier.maxHarmonics, floor(nH/2));
    end
end

y = resp;
if detrendY
    y = y - mean(y, 'omitnan');
end

if useHann
    w = hann(nH);
    y = y .* w;
end

Y = fft(y);
DC = Y(1);
pos = Y(2:floor(nH/2)+1);

harmonics = (1:numel(pos)).';

amp = abs(pos) ./ nH;
phase_deg = rad2deg(angle(pos));

powerAll = sum(abs(pos).^2);
if powerAll > 0
    powerFrac = (abs(pos).^2) ./ powerAll;
else
    powerFrac = nan(size(pos));
end

K = min(maxHarm, numel(harmonics));
harmonics  = harmonics(1:K);
amp        = amp(1:K);
phase_deg  = phase_deg(1:K);
powerFrac  = powerFrac(1:K);

FT.hues        = hues;
FT.resp        = resp;
FT.DC          = DC;
FT.harmonics   = harmonics;
FT.amp         = amp;
FT.phase_deg   = phase_deg;
FT.powerFrac   = powerFrac;
FT.useHann     = useHann;
FT.detrend     = detrendY;

end
