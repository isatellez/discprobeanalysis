function H = plot_vss_unit_harmonics(U, hm, P, COL, config, outDir)

H = struct('fig',[], 'ax',[]);

if U.nTrials == 0
    return;
end

makePlots = true;
if isfield(config,'plot') && isfield(config.plot,'makePlots')
    makePlots = logical(config.plot.makePlots);
end

figVis = 'on';
if ~makePlots
    figVis = 'off';
end

pngDpi = 300;
if isfield(config,'plot') && isfield(config.plot,'dpi')
    pngDpi = config.plot.dpi;
end

dateStrStr  = string(U.dateStr);
unitTypeStr = string(U.unitType);
unitIDVal   = U.unitID;
unitLabel   = unitTypeStr + unitIDVal;    % "SU13", "MU269"

if nargin < 2 || isempty(hm)
    hm = compute_hue_means(U);
end

hues = hm.hues(:);
y    = hm.rate_mean(:);

good = isfinite(hues) & isfinite(y);
hues = hues(good);
y    = y(good);

[hues, idx] = sort(hues);
y = y(idx);

nH = numel(hues);
if nH < 4
    return;
end

cols_satur2 = [];
if isfield(COL,'colsSatur') && ~isempty(COL.colsSatur)
    sz = size(COL.colsSatur);          % nSat x nHue x 3
    cols_satur2 = reshape(double(COL.colsSatur), sz(1), sz(2), 3);
end

cHue = [];
if ~isempty(cols_satur2)
    satRow = size(cols_satur2,1);      % highest saturation row
    c = squeeze(cols_satur2(satRow, hues, :));
    c = double(c);
    if max(c(:)) > 1
        c = c/255;
    end
    cHue = c;
end

F = fit_first_second_harmonics(y);

peakClass = "";
if nargin >= 3 && ~isempty(P) && isfield(P,'class')
    peakClass = string(P.class);
end

ttl = sprintf('%s%d %s %s', char(unitTypeStr), unitIDVal, ...
    char(peakClass), char(dateStrStr));

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

unitFigDir = fullfile(outDir, char(unitLabel));
if ~exist(unitFigDir, 'dir')
    mkdir(unitFigDir);
end

fig = figure('Visible', figVis, 'Color', 'w', 'Position', [100 100 900 600]);
H.fig = fig;

ax = axes('Parent', fig);
H.ax = ax;
hold(ax, 'on');

b = bar(ax, 1:nH, y);
if ~isempty(cHue)
    b.FaceColor = 'flat';
    b.CData     = cHue;
end

xlabel(ax, 'Hue index');
ylabel(ax, 'Mean rate (spikes/s)');
title(ax, ttl, 'Interpreter','none');
grid(ax, 'on');
xlim(ax, [0.5, nH+0.5]);

hx = linspace(1, nH, 400);

y1 = F.a0 + F.amp1 * cos(2*pi*(hx-1)/nH + deg2rad(F.ph1));
y2 = F.a0 + F.amp2 * cos(4*pi*(hx-1)/nH + deg2rad(F.ph2));

if peakClass == "bimodal"
    plot(ax, hx, y1, 'k--', 'LineWidth', 1.5);
    plot(ax, hx, y2, 'k:',  'LineWidth', 1.5);
    legend(ax, {'data','1st harmonic','2nd harmonic'}, 'Location','best');
else
    plot(ax, hx, y1, 'k--', 'LineWidth', 1.5);
    legend(ax, {'data','1st harmonic'}, 'Location','best');
end

txtBox = sprintf(['R^2_1 = %.2f\nR^2_{1+2} = %.2f\n' ...
                  'P_1 = %.2f\nP_2 = %.2f'], ...
                 F.R2_1, F.R2_12, F.pow1, F.pow2);

annotation(fig, 'textbox', [0.65 0.65 0.3 0.25], ...
    'String', txtBox, ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top', ...
    'EdgeColor','none', ...
    'Interpreter','tex');

fileTag = sprintf('%s_%s_%d', char(dateStrStr), char(unitTypeStr), unitIDVal);
pngName = fullfile(unitFigDir, sprintf('VSS_Harmonics_%s.png', fileTag));
figName = fullfile(unitFigDir, sprintf('VSS_Harmonics_%s.fig', fileTag));

exportgraphics(fig, pngName, 'Resolution', pngDpi);
savefig(fig, figName);

if ~makePlots
    close(fig);
end

end


function F = fit_first_second_harmonics(y)
y = y(:);
n = numel(y);
h = (1:n).';

X1  = [ones(n,1), cos(2*pi*(h-1)/n),   sin(2*pi*(h-1)/n)];
X12 = [ones(n,1), cos(2*pi*(h-1)/n),   sin(2*pi*(h-1)/n), ...
                 cos(4*pi*(h-1)/n),   sin(4*pi*(h-1)/n)];

beta1  = X1  \ y;
beta12 = X12 \ y;

a0   = beta12(1);
B1   = beta12(2);
C1   = beta12(3);
B2   = beta12(4);
C2   = beta12(5);

amp1 = sqrt(B1.^2 + C1.^2);
amp2 = sqrt(B2.^2 + C2.^2);

ph1  = atan2(-C1, B1);
ph2  = atan2(-C2, B2);
ph1  = rad2deg(ph1);
ph2  = rad2deg(ph2);

yhat1  = X1  * beta1;
yhat12 = X12 * beta12;

mu    = mean(y);
SStot = sum((y - mu).^2);

if SStot > 0
    R2_1  = 1 - sum((y - yhat1 ).^2) / SStot;
    R2_12 = 1 - sum((y - yhat12).^2) / SStot;
else
    R2_1  = NaN;
    R2_12 = NaN;
end

y1_only = a0 + X1(:,2:3) * beta12(2:3);
y2_only = a0 + [cos(4*pi*(h-1)/n), sin(4*pi*(h-1)/n)] * beta12(4:5);

if SStot > 0
    pow1 = sum((y1_only - mu).^2) / SStot;
    pow2 = sum((y2_only - mu).^2) / SStot;
else
    pow1 = NaN;
    pow2 = NaN;
end

F = struct();
F.a0    = a0;
F.amp1  = amp1;
F.amp2  = amp2;
F.ph1   = ph1;
F.ph2   = ph2;
F.R2_1  = R2_1;
F.R2_12 = R2_12;
F.pow1  = pow1;
F.pow2  = pow2;

end
