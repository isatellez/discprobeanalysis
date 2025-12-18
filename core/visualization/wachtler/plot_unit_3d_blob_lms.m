function H = plot_unit_3d_blob_lms(U, tuning3D, config)
% One 3D response blob per unit in [L-M, S, Lum] space.

H = struct('fig',[], 'ax',[]);

LM  = tuning3D.LM(:);
S   = tuning3D.S(:);
Lum = tuning3D.Lum(:);
R   = tuning3D.rate(:);

ok = isfinite(LM) & isfinite(S) & isfinite(Lum) & isfinite(R);
LM  = LM(ok);
S   = S(ok);
Lum = Lum(ok);
R   = R(ok);

if isempty(LM)
    return;
end

% normalize responses 0–1
R = R - min(R);
if max(R) > 0
    R = R ./ max(R);
end

% grid limits with padding
pad = 0.05;
lm_grid  = linspace(min(LM)-pad,  max(LM)+pad, 20);
s_grid   = linspace(min(S)-pad,   max(S)+pad,  20);
lum_grid = linspace(min(Lum)-pad, max(Lum)+pad,20);

% IMPORTANT: use meshgrid (what isonormals/interp3 expect)
[LMg, Sg, Lumg] = meshgrid(lm_grid, s_grid, lum_grid);

% interpolate responses onto grid
F  = scatteredInterpolant(LM, S, Lum, R, 'natural','none');
Rg = F(LMg, Sg, Lumg);

H.fig = figure('Color','w');
H.ax  = axes('Parent',H.fig);
ax    = H.ax;
hold(ax,'on');

levels = [0.3 0.6 0.9];

for k = 1:numel(levels)
    lvl = levels(k);
    try
        pch = patch(isosurface(LMg, Sg, Lumg, Rg, lvl));
    catch
        continue;
    end
    % skip isonormals to avoid interp3 headaches; surfaces still look fine
    pch.FaceColor = [0.7 0.7 0.7];
    pch.EdgeColor = 'none';
    pch.FaceAlpha = 0.25 + 0.15*k;
end

% overlay sample points
scatter3(ax, LM, S, Lum, 40, R, 'filled', 'MarkerEdgeColor','k');

xlabel(ax,'L - M');
ylabel(ax,'S');
zlabel(ax,'Lum (L + M)');
grid(ax,'on');
axis(ax,'equal');
view(ax, 35, 20);
box(ax,'on');

ttl = sprintf('%s %s %d', string(U.dateStr), U.unitType, U.unitID);
title(ax, ttl);

end
