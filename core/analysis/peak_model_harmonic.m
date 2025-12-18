function P = peak_model_harmonic(U, config)
% Peak-shape classifier using explicit harmonic models:
% flat vs 1st harmonic vs 1st+2nd harmonic

hm   = compute_hue_means(U);
hues = hm.hues;
y    = hm.rate_mean;

P = struct();
P.hues = hues;
P.resp = y;

if nargin < 2
    config = struct();
end

alpha = 0.05;
if isfield(config,'tuning') && isfield(config.tuning,'alphaHarmonic') ...
        && ~isempty(config.tuning.alphaHarmonic)
    alpha = config.tuning.alphaHarmonic;
end

% empty or all NaN → call it flat / undefined stats
if isempty(y) || all(~isfinite(y))
    P.class      = "flat";
    P.p_tuned    = NaN;
    P.p_bimodal  = NaN;
    P.r2_flat    = NaN;
    P.r2_1       = NaN;
    P.r2_12      = NaN;
    P.SSE0       = NaN;
    P.SSE1       = NaN;
    P.SSE2       = NaN;
    P.amp1       = NaN;
    P.amp2       = NaN;
    P.amp_ratio  = NaN;
    P.R2         = NaN;
    return;
end

y = y(:);
n = numel(y);

% not enough directions to do anything sensible
if n < 4
    P.class      = "flat";
    P.p_tuned    = NaN;
    P.p_bimodal  = NaN;
    P.r2_flat    = NaN;
    P.r2_1       = NaN;
    P.r2_12      = NaN;
    P.SSE0       = NaN;
    P.SSE1       = NaN;
    P.SSE2       = NaN;
    P.amp1       = NaN;
    P.amp2       = NaN;
    P.amp_ratio  = NaN;
    P.R2         = NaN;
    return;
end

% drop NaNs if any
isGood = isfinite(y);
if ~all(isGood)
    y = y(isGood);
    n = numel(y);
end

% equally spaced directions around the circle
theta = 2*pi*(0:n-1)'/n;

% design matrices
X0 = ones(n,1);                                % flat
X1 = [ones(n,1), cos(theta), sin(theta)];      % 1st harmonic
X2 = [X1, cos(2*theta), sin(2*theta)];         % 1st + 2nd

% fits
beta0 = X0 \ y;
beta1 = X1 \ y;
beta2 = X2 \ y;

yhat0 = X0*beta0;
yhat1 = X1*beta1;
yhat2 = X2*beta2;

res0 = y - yhat0;
res1 = y - yhat1;
res2 = y - yhat2;

SSE0 = sum(res0.^2);
SSE1 = sum(res1.^2);
SSE2 = sum(res2.^2);

yMean = mean(y);
SST   = sum((y - yMean).^2);

if SST > 0
    r2_flat = 1 - SSE0 / SST;
    r2_1    = 1 - SSE1 / SST;
    r2_12   = 1 - SSE2 / SST;
else
    r2_flat = NaN;
    r2_1    = NaN;
    r2_12   = NaN;
end

% model dimensions
p0 = 1;   % flat: intercept
p1 = 3;   % flat + cos + sin
p2 = 5;   % flat + cos + sin + cos2 + sin2

df0 = n - p0;
df1 = n - p1;
df2 = n - p2;

% F-test: flat vs 1st harmonic
if df1 > 0 && (SSE1 < SSE0)
    num01 = (SSE0 - SSE1) / (p1 - p0);
    den01 = SSE1 / df1;
    F01   = num01 / den01;
    p_tuned = 1 - fcdf(F01, p1 - p0, df1);
else
    p_tuned = 1;
end

% F-test: 1st vs 1st+2nd harmonic
if df2 > 0 && (SSE2 < SSE1)
    num12 = (SSE1 - SSE2) / (p2 - p1);
    den12 = SSE2 / df2;
    F12   = num12 / den12;
    p_bimodal = 1 - fcdf(F12, p2 - p1, df2);
else
    p_bimodal = 1;
end

ratioThresh = 0.4;


% quick FFT amplitudes (for compatibility / extra descriptors)
Yfft = fft(y - min(y));
pos  = Yfft(2:floor(n/2)+1);
amp1 = 0;
amp2 = 0;
if ~isempty(pos)
    amp1 = 2*abs(pos(1)) / n;
    if numel(pos) >= 2
        amp2 = 2*abs(pos(2)) / n;
    end
end
if amp1 > 0
    amp_ratio = amp2 / amp1;
else
    amp_ratio = NaN;
end


ratioThresh = 0.4;

% classification: first decide tuned vs flat from p_tuned,
% then upgrade unimodal → bimodal based on 2nd harmonic
if p_tuned >= alpha
    % not significantly tuned overall → flat
    cls = "flat";
else
    % significantly tuned overall → at least unimodal
    if (p_bimodal < alpha) %|| (amp_ratio >= ratioThresh)
        cls = "bimodal";
    else
        cls = "unimodal";
    end
end


P.class      = cls;
P.p_tuned    = p_tuned;
P.p_bimodal  = p_bimodal;
P.r2_flat    = r2_flat;
P.r2_1       = r2_1;
P.r2_12      = r2_12;
P.SSE0       = SSE0;
P.SSE1       = SSE1;
P.SSE2       = SSE2;
P.amp1       = amp1;
P.amp2       = amp2;
P.amp_ratio  = amp_ratio;
P.R2         = r2_12;   % full model R²
