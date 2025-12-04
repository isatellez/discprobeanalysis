function RY = rayleigh_test(U, config, COL)
% Rayleigh test on hue tuning (early window)

hm = compute_hue_means(U);
hues = hm.hues;
resp = hm.rate_mean;

if all(~isfinite(resp)) || all(resp <= 0 | isnan(resp))
    RY = struct();
    RY.hues      = hues;
    RY.theta_deg = nan(size(hues));
    RY.rate_mean = resp;
    RY.R         = NaN;
    RY.p         = NaN;
    RY.z         = NaN;
    RY.mu_deg    = NaN;
    RY.n         = numel(hues);
    return;
end

resp(resp < 0) = 0;

theta_deg = spatial_transforms.hue_to_angle(hues, COL, config);
theta = deg2rad(theta_deg);

w = resp(:);
C = sum(w .* cos(theta));
S = sum(w .* sin(theta));
Rvec = sqrt(C.^2 + S.^2);
R = Rvec ./ sum(w);

n = numel(hues);
z = n * R.^2;
p = exp(-z);

mu = atan2(S, C);
mu_deg = mod(rad2deg(mu), 360);

RY.hues      = hues;
RY.theta_deg = theta_deg;
RY.rate_mean = resp;
RY.R         = R;
RY.z         = z;
RY.p         = p;
RY.mu_deg    = mu_deg;
RY.n         = n;

end
