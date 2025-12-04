function COL = load_color_metadata(config)

% load ProbeCols + ProbeColIDs
colsFile = fullfile(config.paths.code, 'DiscProbeColsUpdated.mat');
if ~isfile(colsFile)
    error('Missing DiscProbeColsUpdated.mat at %s', colsFile);
end

S = load(colsFile);
COL.probeCols   = double(S.ProbeCols);
COL.probeIDs    = double(S.ProbeColIDs);

COL.nHue = max(COL.probeIDs(:,1));

% reshape into [sat × hue × 3] for saturations
COL.colsSatur = reshape(COL.probeCols(1:48,:), 3, COL.nHue, 3);
COL.colsSatur = sanitize_colors(COL.colsSatur);

% elevation block if present
if size(COL.probeCols,1) >= 138
    COL.colsElev = reshape(COL.probeCols(49:138,:), 6, 15, 3);
else
    COL.colsElev = [];
end

% load DKL CSV
dklPaths = {
    fullfile(config.paths.code,'dkl.csv')
    fullfile(config.paths.base,'dkl.csv')
    'dkl.csv'
};

DKL = [];
for i = 1:numel(dklPaths)
    f = dklPaths{i};
    if isfile(f)
        try
            A = readmatrix(f);
            if size(A,2) >= 3
                DKL = A(:,1:3);
                break;
            end
        catch
        end
    end
end

if isempty(DKL)
    error('Could not load dkl.csv with 3 numeric columns.');
end

% reshape: [sat × hue × axis]
COL.dkl = reshape(DKL(1:48,:), 3, COL.nHue, 3);

end

% small helper
function A = sanitize_colors(A)
    A = double(A);
    A(~isfinite(A)) = 0;
    A = min(max(A,0),255);
end
