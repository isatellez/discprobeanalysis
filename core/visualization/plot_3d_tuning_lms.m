function ax = plot_3d_tuning_lms(tuning3D, varargin)
% Plot 3D tuning contour in [L-M, S, Lum] space,
% similar to Wachtler Fig 1C/D.

p = inputParser;
addParameter(p,'Levels',[0.3 0.6 0.9]);  % iso-response levels (fraction of max)
addParameter(p,'Axes',[]);
addParameter(p,'Title','');
parse(p,varargin{:});

levelsFrac = p.Results.Levels;
ax         = p.Results.Axes;
ttl        = p.Results.Title;

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
    warning('No valid data for 3D tuning.');
    return;
end

% normalize so max response = 1
R = R - min(R);
if max(R) > 0
    R = R ./ max(R);
end

pad = 0.05;
lm_grid  = linspace(min(LM)-pad,  max(LM)+pad, 30);
s_grid   = linspace(min(S)-pad,   max(S)+pad,  30);
lum_grid = linspace(min(Lum)-pad, max(Lum)+pad,30);

[LMg, Sg, Lumg] = ndgrid(lm_grid, s_grid, lum_grid);

F  = scatteredInterpolant(LM, S, Lum, R, 'natural','none');
Rg = F(LMg, Sg, Lumg);

if isempty(ax) || ~isvalid(ax)
    figure('Color','w');
    ax = axes;
end
axes(ax);
hold(ax,'on');

levels = levelsFrac;

for k = 1:numel(levels)
    lvl = levels(k);
    try
        pch = patch(isosurface(LMg, Sg, Lumg, Rg, lvl));
    catch
        continue;
    end
    isonormals(LMg, Sg, Lumg, Rg, pch);
    pch.FaceColor = [0.7 0.7 0.7];
    pch.EdgeColor = 'none';
    pch.FaceAlpha = 0.25 + 0.15*k;  % more opaque at higher response
end

scatter3(ax, LM, S, Lum, 50, R, 'filled', 'MarkerEdgeColor','k');

xlabel(ax,'L - M');
ylabel(ax,'S');
zlabel(ax,'Lum (L + M)');

grid(ax,'on');
axis(ax,'equal');
view(ax, 30, 20);
box(ax,'on');

if ~isempty(ttl)
    title(ax, ttl);
end

end
