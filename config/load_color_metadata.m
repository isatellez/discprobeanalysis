function COL = load_color_metadata(config)

% load ProbeCols + ProbeColIDs
colsFile = fullfile(config.paths.code, 'DiscProbeColsUpdated.mat');
% this file has hue, saturation, elevation, dkl coords, rgb values
% for the 144 stimuli used in the task

if ~isfile(colsFile)
    error('Missing DiscProbeColsUpdated.mat at %s', colsFile);
end

S = load(colsFile);

COL = struct();

COL.probeCols = sanitize_colors(S.ProbeCols);   % nProbes × 3 RGB
COL.probeIDs  = double(S.ProbeColIDs);          % nProbes × 6: [hue sat elev, dkl]

rgbAll = COL.probeCols;
ids    = COL.probeIDs;

hueID  = ids(:,1); % hue id (should be integer-like)
satID  = ids(:,2); % 0.2, 0.33, 0.5, 0.66, 0.8, 1 or similar
elevID = ids(:,3); % elevation values (e.g. -0.8 ... 0.8)

% make sure hue IDs are integer-ish
hueID_int = round(hueID);
if any(abs(hueID_int - hueID) > 1e-6)
    warning('load_color_metadata: hueID values are not close to integers.');
end
hueID = hueID_int;

COL.nHue = max(hueID(:));

% saturation / elevation levels are value-based, not 1..N indexes
satVals  = unique(satID(~isnan(satID)));
elevVals = unique(elevID(~isnan(elevID)));

COL.nSat      = numel(satVals);
COL.nElev     = numel(elevVals);
COL.satVals   = satVals;
COL.elevVals  = elevVals;

nHue  = COL.nHue;
nSat  = COL.nSat;
nElev = COL.nElev;

% make sure RGB are clean and preallocate lookups
COL.colsSatur = nan(nSat,  nHue, 3);
COL.colsElev  = nan(nElev, nHue, 3);

for k = 1:numel(hueID)
    h   = hueID(k);
    sVal = satID(k);
    eVal = elevID(k);

    if ~isfinite(h) || h < 1 || h > nHue
        continue;
    end

    % map satID/elevID values to index in satVals/elevVals
    sIdx = find(satVals == sVal, 1);
    eIdx = find(elevVals == eVal, 1);

    if ~isempty(sIdx)
        COL.colsSatur(sIdx, h, :) = rgbAll(k,:);
    end
    if ~isempty(eIdx)
        COL.colsElev(eIdx, h, :) = rgbAll(k,:);
    end
end

% load DKL CSV 
dklFile = fullfile(config.paths.code, 'dkl.csv');

if ~isfile(dklFile)
    error('Could not find dkl.csv at %s', dklFile);
end

A = readmatrix(dklFile);
if size(A,2) < 3
    error('dkl.csv must have at least 3 numeric columns, found %d.', size(A,2));
end

DKL = A(:,1:3);

% sanity check row counts
nProbes = size(COL.probeCols,1);
if size(DKL,1) ~= nProbes
    error('dkl.csv row count (%d) does not match ProbeCols (%d).', ...
          size(DKL,1), nProbes);
end

% row-wise DKL for master-row indexing
COL.dklRows = DKL;   % <-- this is what TrialIndex + index_to_dkl will use

% build DKL lookup by saturation × hue (for DKL-plane utilities)
COL.dkl = nan(nSat, nHue, 3);
for k = 1:numel(hueID)
    h   = hueID(k);
    sVal = satID(k);

    if ~isfinite(h) || h < 1 || h > nHue
        continue;
    end

    sIdx = find(satVals == sVal, 1);
    if ~isempty(sIdx)
        COL.dkl(sIdx, h, :) = DKL(k,:);
    end
end

% load LMS CSV 
lmsFile = fullfile(config.paths.code, 'lms.csv');

LMS = [];

if isfile(lmsFile)
    try
        A = readmatrix(lmsFile);
        if size(A,2) < 3
            warning('lms.csv must have at least 3 numeric columns, found %d. LMS will be unavailable.', size(A,2));
        else
            LMS = A(:,1:3);
        end
    catch ME
        warning('Failed to read lms.csv at %s: %s. LMS will be unavailable.', ...
                lmsFile, ME.message);
    end
else
    warning('Could not find lms.csv at %s. LMS will be unavailable.', lmsFile);
end

if ~isempty(LMS)
    if size(LMS,1) ~= nProbes
        warning('lms.csv row count (%d) does not match ProbeCols (%d).', ...
                size(LMS,1), nProbes);
    end

    COL.lms = LMS;

    L     = LMS(:,1);
    M     = LMS(:,2);
    Scone = LMS(:,3);

    denom = L + M;
    denom(denom == 0) = NaN;

    COL.lmsCC = [L./denom, M./denom, Scone./denom]; %normalized lms 
end

end

function A = sanitize_colors(A) %to clip any values that go below or above 0-255
    A = double(A);
    A(~isfinite(A)) = 0;
    A = min(max(A,0),255);
end
