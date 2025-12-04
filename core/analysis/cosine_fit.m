function COS = cosine_fit(U, config)
% Cosine fit to hue tuning (first harmonic model)

hm = compute_hue_means(U);
hues = hm.hues;
y    = hm.rate_mean;

COS = struct();
COS.hues      = hues;
COS.theta_deg = nan(size(hues));
COS.y         = y;
COS.y_hat     = nan(size(y));
COS.beta      = nan(3,1);
COS.amp       = NaN;
COS.pref_deg  = NaN;
COS.R2        = NaN;

if all(~isfinite(y)) || all(y == 0 | isnan(y))
    return;
end

theta_deg = spatial_transforms.hue_to_angle(hues, [], config);
theta = theta_deg(:);
y = y(:);

X = [ones(size(theta)) cosd(theta) sind(theta)];
beta = X \ y;

y_hat = X * beta;

amp = hypot(beta(2), beta(3));
phi = atan2(-beta(3), beta(2));
pref_deg = mod(rad2deg(phi), 360);

SS_tot = sum((y - mean(y, 'omitnan')).^2, 'omitnan');
SS_res = sum((y - y_hat).^2, 'omitnan');
if SS_tot > 0
    R2 = 1 - SS_res / SS_tot;
else
    R2 = NaN;
end

COS.theta_deg = theta_deg;
COS.y_hat     = y_hat;
COS.beta      = beta;
COS.amp       = amp;
COS.pref_deg  = pref_deg;
COS.R2        = R2;

end
