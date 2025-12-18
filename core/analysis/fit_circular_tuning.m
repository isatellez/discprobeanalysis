function fit = fit_circular_tuning(thetaDeg, resp, opts)
% Fit Wachtler-style circular normal tuning curve
% fc = A0 + A * exp((cos(phi - phi0) - 1) / sigma^2)
% Also returns chi-square goodness-of-fit and p-value.

    theta = deg2rad(thetaDeg(:));
    y     = resp(:);

    % remove NaNs
    valid = isfinite(theta) & isfinite(y);
    theta = theta(valid);
    y     = y(valid);

    fit = struct('A0',NaN,'A',NaN,'phi0_deg',NaN,'sigma',NaN, ...
                 'fwhm_deg',NaN,'R2',NaN,'chi2',NaN,'pChi2',NaN,'df',NaN);

    if numel(y) < 4
        return;
    end

    if nargin < 3
        opts = struct();
    end

    if ~isfield(opts,'A0_init')
        opts.A0_init = min(y);
    end
    if ~isfield(opts,'A_init')
        opts.A_init = max(y) - opts.A0_init;
    end

    [~, idxMax] = max(y);
    if ~isfield(opts,'phi0_init')
        opts.phi0_init = theta(idxMax);
    end
    if ~isfield(opts,'sigma_init')
        opts.sigma_init = 1;
    end

    p0 = [opts.A0_init, opts.A_init, opts.phi0_init, opts.sigma_init];

    % objective: sum of squared errors
    obj = @(p) sse_circ_model(p, theta, y);

    opts_fmin = optimset('Display','off');
    p_hat = fminsearch(obj, p0, opts_fmin);

    A0   = p_hat(1);
    A    = p_hat(2);
    phi0_rad = p_hat(3);
    phi0 = phi0_rad;
    sig  = max(p_hat(4), 1e-3);
    % A0   = p_hat(1);
    % A    = p_hat(2);
    % phi0 = p_hat(3);
    % sig  = max(p_hat(4), 1e-3);

    yhat = circ_model(p_hat, theta);

    % R^2
    ss_res = sum((y - yhat).^2);
    ymean  = mean(y);
    ss_tot = sum((y - ymean).^2);
    if ss_tot <= eps
        R2 = NaN;
    else
        R2 = 1 - ss_res / ss_tot;
    end

    % FWHM in deg
    cos_dphi = 1 - sig.^2 * log(2);
    cos_dphi = min(max(cos_dphi, -1), 1);
    dphi_rad = acos(cos_dphi);
    fwhm_deg = rad2deg(2 * dphi_rad);

    % chi-square goodness-of-fit
    % Use Poisson-like variance: var ~ mean rate
    varY = max(yhat, 1e-3);
    chi2 = sum((y - yhat).^2 ./ varY);

    k  = 4;                 % A0, A, phi0, sigma
    N  = numel(y);
    df = max(N - k, 1);     % for 8 dirs -> df = 4

    pChi2 = 1 - chi2cdf(chi2, df);

    fit.A0       = A0;
    fit.A        = A;
    fit.phi0_deg = mod(rad2deg(phi0), 360);
    fit.sigma    = sig;
    fit.fwhm_deg = fwhm_deg;
    fit.R2       = R2;
    fit.chi2     = chi2;
    fit.pChi2    = pChi2;
    fit.df       = df;
end

function yhat = circ_model(p, theta)
% circular normal model

    A0   = p(1);
    A    = p(2);
    phi0 = p(3);
    sig  = max(p(4), 1e-3);

    yhat = A0 + A * exp((cos(theta - phi0) - 1) ./ (sig.^2));
end

function sse = sse_circ_model(p, theta, y)
% sum of squared errors between data and model

    yhat = circ_model(p, theta);
    sse  = sum((y - yhat).^2);
end
