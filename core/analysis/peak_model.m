function P = peak_model(U, config)
% simple peak-shape classifier based on Fourier harmonics
% labels: 'flat', 'unimodal', 'bimodal'

hm = compute_hue_means(U);
hues = hm.hues;
y    = hm.rate_mean;

P = struct();
P.hues = hues;
P.resp = y;

if all(~isfinite(y)) || all(y <= 0 | isnan(y))
    P.class     = "flat";
    P.amp1      = NaN;
    P.amp2      = NaN;
    P.amp_ratio = NaN;
    P.p_tuned   = NaN;
    P.p_bimodal = NaN;
    return;
end

y = y(:);
y = y - min(y);
if max(y) > 0
    y = y ./ max(y);
end

nH = numel(y);
if nH < 4
    P.class     = "flat";
    P.amp1      = NaN;
    P.amp2      = NaN;
    P.amp_ratio = NaN;
    P.p_tuned   = NaN;
    P.p_bimodal = NaN;
    return;
end

Y = fft(y);
pos = Y(2:floor(nH/2)+1);

amp1 = 2*abs(pos(1)) ./ nH;
if numel(pos) >= 2
    amp2 = 2*abs(pos(2)) ./ nH;
else
    amp2 = 0;
end

if amp1 > 0
    amp_ratio = amp2 ./ amp1;
else
    amp_ratio = NaN;
end

if amp1 < 0.1
    cls = "flat";
elseif amp_ratio > 0.5
    cls = "bimodal";
else
    cls = "unimodal";
end

P.class     = cls;
P.amp1      = amp1;
P.amp2      = amp2;
P.amp_ratio = amp_ratio;

P.p_tuned   = NaN;   % we can hook in permutation/Fourier p's later
P.p_bimodal = NaN;

end
