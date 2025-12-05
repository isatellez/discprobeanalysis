function plot_dkl_suite(U, stats, COL, config, outDir)

if U.nTrials == 0
    return;
end

if ~isfield(COL,'dklSatur') || isempty(COL.dklSatur)
    % need DKL per sat×hue from color metadata
    return;
end

dkl_satur2  = COL.dklSatur;   % nSat×nHue×3, legacy style
cols_satur2 = COL.colsSatur;  % nSat×nHue×3 RGB

nSat_meta = size(dkl_satur2,1);
nHue_meta = size(dkl_satur2,2);

% basic IDs from trials
hueID = U.trials.hueID;
satID = U.trials.satID;

saturIDs = unique(satID(~isnan(satID)));
saturIDs = saturIDs(:).';
nS       = numel(saturIDs);

% map trial sats into 1..nSat indices
[~, satIdx_u] = ismember(satID, saturIDs);
hueIdx_u      = hueID;

% early window in Hz
win = U.winEarly;
dur = diff(win);

nTr = numel(U.spk);
spk_early_hz = nan(nTr,1);
for tt = 1:nTr
    spks = U.spk{tt};
    if isempty(spks)
        spk_early_hz(tt) = 0;
    else
        spks = spks(spks >= win(1) & spks < win(2));
        spk_early_hz(tt) = numel(spks) / dur;
    end
end

% mean response per (sat,hue) like resp_by_sh in the script
R = NaN(nS, nHue_meta);
for si = 1:nS
    for h = 1:nHue_meta
        sel = (satIdx_u == si) & (hueIdx_u == h);
        R(si,h) = mean(spk_early_hz(sel), 'omitnan');
    end
end

% if metadata has more sats than data, trim; if fewer, cap
if nS ~= nSat_meta
    m = min(nS, nSat_meta);
    R           = R(1:m,:);
    saturIDs    = saturIDs(1:m);
    dkl_satur2  = dkl_satur2(1:m,:,:);
    cols_satur2 = cols_satur2(1:m,:,:);
    nS          = m;
    nSat_meta   = m;
end

nHue = nHue_meta;

% detect "L+M" axis as lowest variance across hue
D = reshape(dkl_satur2(:,1:nHue,:), [], 3);   % (nSat*nHue)×3
score = std(D,0,1) + 0.1*mean(abs(D),1);
[~, idxLum] = min(score);
idxCh = setdiff(1:3, idxLum);
ixRG  = idxCh(1);
ixS   = idxCh(2);

% isoluminant plane coords per (sat,hue)
XY = reshape(dkl_satur2(:,:,[ixRG ixS]), [], 2);   % (nSat*nHue)×2

% colors for each bar
C = double(reshape(cols_satur2, [], 3));
if isempty(C) || size(C,2) ~= 3
    C = repmat([0.5 0.5 0.5], size(XY,1), 1);
end
if max(C(:)) > 1
    C = C / 255;
end
C(~isfinite(C)) = 0;
C = min(max(C,0),1);

% flatten response
r = reshape(R, [], 1);

% mask out bad points
bad = any(~isfinite(XY),2) | ~isfinite(r);
XY(bad,:) = [];
C(bad,:)  = [];
r(bad)    = [];

if isempty(XY) || isempty(r)
    return;
end

% mean per hue in plane, for the polar-ish summary
XY_all = reshape(dkl_satur2(:,:,[ixRG ixS]), [], 2);
XY_all = reshape(XY_all, [nSat_meta*nHue, 2]);
R_all  = reshape(R, [], 1);
good   = all(isfinite(XY_all),2) & isfinite(R_all);
XY_all = XY_all(good,:);
R_all  = R_all(good);

if isempty(XY_all)
    XYmean = nan(nHue,2);
else
    XYmean = nan(nHue,2);
    for h = 1:nHue
        rows = ((0:nSat_meta-1)*nHue) + h;
        rows = rows(rows <= size(XY_all,1));
        XYmean(h,:) = nanmean(XY_all(rows,:),1);
    end
end

th     = mod(atan2(XYmean(:,2), XYmean(:,1)), 2*pi);
r_hue  = max(0, nanmean(R,1).');
[thS, ord] = sort(th);
rS         = r_hue(ord);

% plotting flags
makePlots = true;
if isfield(config,'plot') && isfield(config.plot,'makePlots')
    makePlots = logical(config.plot.makePlots);
end
figVis = 'on';
if ~makePlots
    figVis = 'off';
end

unitDir = fullfile(outDir, sprintf('%s_%d', U.unitType, U.unitID));
if ~exist(unitDir,'dir'), mkdir(unitDir); end

f = figure('Color','w','Visible',figVis,'Name','DKL suite');
f.Units    = 'normalized';
f.Position = [0.1 0.1 0.75 0.7];

tlo = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');

% 3D bars in DKL chromatic plane
ax3d = nexttile(tlo,1);
hold(ax3d,'on'); grid(ax3d,'on'); axis(ax3d,'equal');

% scale radius to something reasonable
rPos = max(r,0);
maxR = max(rPos);
if maxR <= 0
    rPos(:) = 0;
end

% bar "footprint" radius relative to max response
baseRadius = 0.06;
for k = 1:numel(rPos)
    if rPos(k) <= 0, continue; end
    x0 = XY(k,1);
    y0 = XY(k,2);
    z0 = 0;
    h  = rPos(k);

    w  = baseRadius;
    xv = [x0-w x0+w x0+w x0-w];
    yv = [y0-w y0-w y0+w y0+w];
    zv = [z0   z0   z0   z0];
    zt = zv + h;

    fc = C(k,:);
    patch(ax3d, xv, yv, zt, fc, 'FaceAlpha',0.9, 'EdgeColor','none');
    patch(ax3d, xv, yv, zv, fc*0.6, 'FaceAlpha',0.7, 'EdgeColor','none');
end

xlabel(ax3d,'DKL chromatic axis 1');
ylabel(ax3d,'DKL chromatic axis 2');
zlabel(ax3d,'Rate (Hz)');
title(ax3d,'Mean response per (sat,hue) in DKL plane');

view(ax3d, [40 25]);

% polar-ish summary of hue tuning in that plane
ax2d = nexttile(tlo,2);
hold(ax2d,'on'); grid(ax2d,'on');

if all(isfinite(thS)) && any(rS>0)
    polaraxes(ax2d.Parent); delete(ax2d); % ensure polar axes
    ax2d = polaraxes('Parent',f);
    ax2d.Layout.Tile = 2;
    polarplot(ax2d, thS, rS, 'o-','LineWidth',1.5);
    title(ax2d,'Mean rate vs hue angle in DKL plane');
else
    % fallback if angles are bad
    plot(ax2d, 1:nHue, r_hue, 'o-','LineWidth',1.5);
    xlim(ax2d,[0.5 nHue+0.5]);
    xlabel(ax2d,'Hue index');
    ylabel(ax2d,'Mean rate (Hz)');
    title(ax2d,'Mean rate per hue');
end

% title
if isfield(U,'exptName')
    ttl = sprintf('%s | Unit %d (%s) — DKL suite', ...
        string(U.exptName), U.unitID, U.unitType);
else
    ttl = sprintf('%s | Unit %d (%s) — DKL suite', ...
        string(U.dateStr), U.unitID, U.unitType);
end
sgtitle(f, ttl, 'Interpreter','none');

% save
fileTag = sprintf('%s_%s_%d', string(U.dateStr), U.unitType, U.unitID);
pngName = fullfile(unitDir, sprintf('DKL_suite_%s.png', fileTag));
figName = fullfile(unitDir, sprintf('DKL_suite_%s.fig', fileTag));

exportgraphics(f, pngName, 'Resolution', 300);
savefig(f, figName);

if ~makePlots
    close(f);
end

end
