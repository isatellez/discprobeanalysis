function H = plot_unit_color_tuning_wachtler(U, unitTuning, COL, config)
% Wachtler-style figure: central polar tuning + surrounding PSTHs.

H = struct('fig',[], 'axPolar',[], 'axPSTH',[]);

thetaDeg = unitTuning.thetaDeg(:);
satIDs   = unitTuning.saturIDs(:).';
meanRate = unitTuning.meanRate;
fit      = [];
if isfield(unitTuning,'fit')
    fit = unitTuning.fit;
end

nH = numel(thetaDeg);
if nH == 0 || isempty(meanRate)
    return;
end

useSats = [];
if isfield(config,'tuning') && isfield(config.tuning,'useSaturIDs') ...
        && ~isempty(config.tuning.useSaturIDs)
    useSats = config.tuning.useSaturIDs(:).';
end

satTol = 1e-3;
if isfield(config,'legacy') && isfield(config.legacy,'satTol')
    satTol = config.legacy.satTol;
end

% choose best saturation for tuning curve
if isempty(satIDs)
    bestSatIdx = 1;
else
    mask = true(size(satIDs));
    if ~isempty(useSats)
        mask = false(size(satIDs));
        for s = useSats
            mask = mask | abs(satIDs - s) < satTol;
        end
    end
    idx = find(mask);
    if isempty(idx)
        idx = 1:numel(satIDs);
    end
    [~, kLocal] = max(satIDs(idx));
    bestSatIdx = idx(kLocal);
end
bestSatIdx = max(1, min(bestSatIdx, size(meanRate,2)));

if ~isempty(satIDs) && bestSatIdx <= numel(satIDs)
    bestSatVal = satIDs(bestSatIdx);
else
    bestSatVal = NaN;
end

resp = meanRate(:, bestSatIdx);
baselineHz = NaN;
if isfield(unitTuning,'baselineHz')
    baselineHz = unitTuning.baselineHz;
end

if all(~isfinite(resp))
    resp(:) = 0;
end

maxR_data = max(resp(:));
if isfinite(baselineHz)
    maxR_data = max(maxR_data, baselineHz);
end

maxR_fit = 0;
if ~isempty(fit)
    phi_tmp = linspace(0,2*pi,360).';
    fitRate_tmp = fit.A0 + fit.A * exp((cos(phi_tmp - deg2rad(fit.phi0_deg)) - 1) / (fit.sigma.^2));
    fitRate_tmp(fitRate_tmp < 0) = 0;
    maxR_fit = max(fitRate_tmp);
end

maxR = max(maxR_data, maxR_fit);
if ~isfinite(maxR) || maxR <= 0
    maxR = 10;
end
maxR = ceil(maxR/10)*10;

thetaRad = deg2rad(thetaDeg);

% ================= central polar axis =================
H.fig = figure('Color','w');
cx = 0.5; cy = 0.5;

polSize = 0.32;
H.axPolar = polaraxes('Parent',H.fig,'Units','normalized', ...
    'Position',[cx-polSize/2 cy-polSize/2 polSize polSize]);

ax = H.axPolar;
hold(ax,'on');
ax.ThetaZeroLocation = 'right';
ax.ThetaDir          = 'counterclockwise';
ax.ThetaTick         = [0 90 180 270];
rlim(ax,[0 maxR]);
rticks(ax,0:20:maxR);

if isfinite(baselineHz) && baselineHz > 0
    phi = linspace(0,2*pi,360).';
    polarplot(ax, phi, baselineHz*ones(size(phi)), 'k:','LineWidth',1);
end

% tuning curve: fit if we have it, else data
if ~isempty(fit)
    phi = linspace(0,2*pi,360).';
    fitRate = fit.A0 + fit.A * exp((cos(phi - deg2rad(fit.phi0_deg)) - 1) / (fit.sigma.^2));
    fitRate(fitRate < 0) = 0;
    polarplot(ax, phi, fitRate, 'k-','LineWidth',3);
else
    respPlot = resp;
    respPlot(~isfinite(respPlot) | respPlot < 0) = 0;
    [thetaSort, ordSort] = sort(thetaRad);
    rSort = respPlot(ordSort);
    polarplot(ax, [thetaSort; thetaSort(1)], [rSort; rSort(1)], ...
        'k-','LineWidth',3);
end

for k = 1:nH
    r = resp(k);
    if ~isfinite(r) || r <= 0
        continue;
    end
    ang = thetaRad(k);
    r0  = max(baselineHz,0);
    if ~isfinite(r0) || r0 < 0
        r0 = 0;
    end
    polarplot(ax, [ang ang], [r0 r], 'k-', 'LineWidth', 1.5);
end

rLbl = maxR * 1.05;
text(ax, deg2rad(90), rLbl, 'S',   'HorizontalAlignment','center', 'VerticalAlignment','bottom');
text(ax, deg2rad(0),  rLbl, 'L-M', 'HorizontalAlignment','left',   'VerticalAlignment','middle');

% ================= title =================
if strcmpi(U.unitType,'SU')
    unitLabel = sprintf('UnitID %d (SU)', U.unitID);
else
    unitLabel = sprintf('MUA ID %d', U.unitID);
end

if isfield(U,'unitIdx')
    unitNumStr = sprintf('UnitNUM %d', U.unitIdx);
elseif isfield(U,'idx')
    unitNumStr = sprintf('UnitNUM %d', U.idx);
else
    unitNumStr = '';
end

if isfield(U,'sessionStamp')
    sessionStamp = string(U.sessionStamp);
elseif isfield(U,'sessionName')
    sessionStamp = sprintf('%s_%s', string(U.dateStr), string(U.sessionName));
elseif isfield(config,'monkey') && ~isempty(config.monkey)
    sessionStamp = sprintf('%s_%s', string(U.dateStr), string(config.monkey));
else
    sessionStamp = string(U.dateStr);
end

prefStr = 'pref=?°';
fwhmStr = 'FWHM=?°';
if ~isempty(fit)
    prefStr = sprintf('pref=%.0f°', fit.phi0_deg);
    fwhmStr = sprintf('FWHM=%.0f°', fit.fwhm_deg);
end

if isfinite(bestSatVal)
    if ~isempty(unitNumStr)
        ttl = sprintf('%s | %s - %s | sat=%.2f | %s | %s', ...
            sessionStamp, unitLabel, unitNumStr, bestSatVal, prefStr, fwhmStr);
    else
        ttl = sprintf('%s | %s | sat=%.2f | %s | %s', ...
            sessionStamp, unitLabel, bestSatVal, prefStr, fwhmStr);
    end
else
    if ~isempty(unitNumStr)
        ttl = sprintf('%s | %s - %s | %s | %s', ...
            sessionStamp, unitLabel, unitNumStr, prefStr, fwhmStr);
    else
        ttl = sprintf('%s | %s | %s | %s', ...
            sessionStamp, unitLabel, prefStr, fwhmStr);
    end
end

title(ax, ttl,'Interpreter','none','FontWeight','bold');
% ===================== PSTHs around the circle =====================

hueTrial = U.trials.hueID;
nTr      = numel(hueTrial);

hasOutliers = isfield(U,'isOutlier') && numel(U.isOutlier)==nTr;

baseMask = isfinite(hueTrial);
if hasOutliers
    baseMask = baseMask & ~U.isOutlier(:);
end

% restrict to isoluminant plane if we have a lum/elev field
if isfield(U.trials,'lumID')
    lumTrial = U.trials.lumID;
    isoVal   = 0;
    lumTol   = 1e-3;
    baseMask = baseMask & abs(lumTrial - isoVal) < lumTol;
elseif isfield(U.trials,'elevID')
    elevTrial = U.trials.elevID;
    isoVal    = 0;
    elevTol   = 1e-3;
    baseMask  = baseMask & abs(elevTrial - isoVal) < elevTol;
end

% unique hues on that plane
huesAvail = unique(hueTrial(baseMask));
huesAvail(isnan(huesAvail)) = [];

if isempty(huesAvail)
    return;
end

% map those hue IDs to angles (same convention as elsewhere)
thetaAvail_deg = spatial_transforms.hue_to_angle(huesAvail, COL, config);
valid = isfinite(thetaAvail_deg);
huesAvail     = huesAvail(valid);
thetaAvail_deg = mod(thetaAvail_deg(valid), 360);

if isempty(huesAvail)
    return;
end

% we WANT 8 PSTHs at 8 target angles, always
nTarget      = 8;
targetTheta  = linspace(0, 360-360/nTarget, nTarget);  % 0,45,...,315

huesShow  = nan(1, nTarget);
thetaShow = targetTheta;   % where to place each PSTH around the circle

for k = 1:nTarget
    % circular distance between available angles and this target
    d = abs(mod(thetaAvail_deg - targetTheta(k) + 180, 360) - 180);
    [~, idxBest] = min(d);
    huesShow(k) = huesAvail(idxBest);   % use nearest hue for this slot
end

nShow = nTarget;

% PSTH binning
binMs   = 10;
tWinMs  = [0 500];
tEdgesMs   = tWinMs(1):binMs:tWinMs(2);
tEdgesSec  = tEdgesMs / 1000;
tCentersMs = tEdgesMs(1:end-1) + binMs/2;
nBins      = numel(tCentersMs);

psthRates = nan(nShow, nBins);

for ii = 1:nShow
    hID = huesShow(ii);
    if isnan(hID)
        continue;
    end

    mask  = baseMask & (hueTrial == hID);
    trIdx = find(mask);
    if isempty(trIdx)
        continue;
    end

    counts = zeros(1,nBins);
    for jj = 1:numel(trIdx)
        sp = U.spk{trIdx(jj)};
        if isempty(sp), continue; end
        counts = counts + histcounts(sp, tEdgesSec);
    end
    psthRates(ii,:) = counts / (numel(trIdx) * (binMs/1000));   % spikes/s
end

% common y-scale
maxPSTH = max(psthRates(:));
if ~isfinite(maxPSTH) || maxPSTH <= 0
    maxPSTH = 1;
end
maxPSTH = ceil(maxPSTH/10)*10;

% geometry: 8 PSTHs around the polar
psthSize = 0.18;
rAxes    = 0.36;

H.axPSTH = gobjects(nShow,1);

for ii = 1:nShow
    th = deg2rad(thetaShow(ii));
    x = cx + rAxes*cos(th) - psthSize/2;
    y = cy + rAxes*sin(th) - psthSize/2;

    H.axPSTH(ii) = axes('Parent',H.fig,'Units','normalized', ...
        'Position',[x y psthSize psthSize]);
    hold(H.axPSTH(ii),'on');

    rates = psthRates(ii,:);
    if any(isfinite(rates))
        bar(H.axPSTH(ii), tCentersMs, rates, 1.0, 'k');
    end
    xlim(H.axPSTH(ii), tWinMs);
    ylim(H.axPSTH(ii), [0 maxPSTH]);
    set(H.axPSTH(ii),'Box','off','TickDir','out','FontSize',7);
    set(H.axPSTH(ii),'XTick',[0 250 500]);
    set(H.axPSTH(ii),'YTick',[0 maxPSTH]);

    if ii ~= 1
        set(H.axPSTH(ii),'XTickLabel',[],'YTickLabel',[]);
    end
end

if ~isempty(H.axPSTH)
    xlabel(H.axPSTH(1),'time [ms]');
    ylabel(H.axPSTH(1),'rate [spikes/s]');
end


end
