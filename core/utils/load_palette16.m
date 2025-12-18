function pal16 = load_palette16(csvPath)
%LOAD_PALETTE16 Load 16 full-saturation RGB colors for the hue ring.
%
%   pal16 = load_palette16(csvPath)
%   - csvPath: path to rgb.csv
%   - pal16:   16x3 matrix of RGB values in [0,1], row 1 is the first hue
%
% Assumes rgb.csv has rows arranged so that every 3rd row is a
% full-saturation color (rows 3,6,9,...).

pal16 = [];

if nargin < 1 || isempty(csvPath)
    warning('load_palette16:NoPath', 'No csvPath provided.');
    return;
end

if ~isfile(csvPath)
    warning('load_palette16:MissingFile', 'rgb.csv not found at %s', csvPath);
    return;
end

A = readmatrix(csvPath);
if isempty(A) || size(A,2) < 3
    warning('load_palette16:BadFile', 'rgb.csv has wrong shape or is empty.');
    return;
end

% every 3rd row is full saturation: 3,6,9,...
idx = 3:3:size(A,1);
if isempty(idx)
    warning('load_palette16:NoFullSatRows', ...
        'No full-saturation rows (3:3:end) found in %s', csvPath);
    return;
end

% pick the first 16 full-saturation hues (or as many as available)
if numel(idx) < 16
    warning('load_palette16:FewHues', ...
        'Only %d full-sat hues found, need 16. Using what is available.', numel(idx));
    idx = idx(1:numel(idx));
else
    idx = idx(1:16);
end

P = double(A(idx,1:3));
P(~isfinite(P)) = 0;

% rgb.csv is effectively 0–1 with occasional 255/256 values
if max(P(:)) > 1
    % assume 0–255 and normalize
    P = P / max(255, max(P(:)));
end

% clip to [0,1] just in case
P = max(0, min(P, 1));

% rotate so the last hue becomes the first (so hue 16 sits at 0°)
if size(P,1) >= 16
    P = circshift(P, 1, 1);   % [16 1 2 ... 15]
end

pal16 = P;
end
