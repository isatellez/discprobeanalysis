function META = build_extended_color_metadata(DATES)
% Scan ExpTrialsDisc files and collect all unique DiscProbe RGBs
% plus which ones are already in DiscProbeColsUpdated and which are new.
%
% This does not change the main analysis. It just writes out a metadata
% .mat and .csv you can use later.

clc;

rootDir = fileparts(mfilename('fullpath'));
addpath(genpath(rootDir));

config = load_config();

if nargin < 1
    DATES = [];
end

monkey = config.monkey;

PATHS = struct();
PATHS.baseDiscProbeLocal = config.paths.base;
PATHS.baseDiscProbesCode = config.paths.code;
PATHS.baseOut            = config.paths.output;

% legacy master from DiscProbeColsUpdated
colsFile = fullfile(config.paths.code, 'DiscProbeColsUpdated.mat');
if ~isfile(colsFile)
    error('Missing DiscProbeColsUpdated.mat at %s', colsFile);
end
S = load(colsFile, 'ProbeCols', 'ProbeColIDs');

legacyRGB = double(S.ProbeCols);    % N0 x 3
legacyIDs = double(S.ProbeColIDs);  % N0 x 6 (hue, sat, elev, ...)

% figure out which dates to scan
if isempty(DATES)
    % auto-detect from any *ExpTrialsDisc.mat under base
    pat = sprintf('%s_*_ExpTrialsDisc.mat', monkey);
    D = dir(fullfile(PATHS.baseDiscProbeLocal, '**', pat));
    if isempty(D)
        error('No ExpTrialsDisc files found under %s', PATHS.baseDiscProbeLocal);
    end
    dateList = strings(numel(D),1);
    for k = 1:numel(D)
        tok = regexp(D(k).name, '\d{6}', 'match', 'once');
        if ~isempty(tok)
            dateList(k) = string(tok);
        end
    end
    dateList = unique(dateList(dateList ~= ""));
else
    dateList = normalize_dates(DATES);
    dateList = string(dateList(:));
end

allRGB      = [];
allDates    = strings(0,1);

fprintf('Scanning %d sessions for DiscprobeColor RGBs...\n', numel(dateList));

for d = 1:numel(dateList)
    dateStr = char(dateList(d));
    try
        expFile = locate_expTrialsDisc(PATHS, monkey, dateStr);
    catch ME
        warning('Skipping %s: %s', dateStr, ME.message);
        continue;
    end

    fprintf('  %s -> %s\n', dateStr, expFile);

    Sexp = load(expFile, 'ExptTrialsDisc');
    if ~isfield(Sexp, 'ExptTrialsDisc')
        warning('  %s: no ExptTrialsDisc in file, skipping.', dateStr);
        continue;
    end

    ExptTrialsDisc = Sexp.ExptTrialsDisc;
    if ~iscell(ExptTrialsDisc) && ~istable(ExptTrialsDisc)
        warning('  %s: ExptTrialsDisc is not a cell/table, skipping.', dateStr);
        continue;
    end

    % pull DiscprobeColor per trial into N x 3 matrix
    try
        if iscell(ExptTrialsDisc)
            trlCols = cell2mat( ...
                cellfun(@(s) double(s.DiscprobeColor(:).'), ...
                        ExptTrialsDisc(:,1), 'UniformOutput', false) );
        else
            % table case: assume a column with structs
            trlCols = cell2mat( ...
                cellfun(@(s) double(s.DiscprobeColor(:).'), ...
                        ExptTrialsDisc.Disc, 'UniformOutput', false) );
        end
    catch ME
        warning('  %s: could not extract DiscprobeColor (%s), skipping.', ...
                dateStr, ME.message);
        continue;
    end

    if isempty(trlCols)
        continue;
    end

    allRGB   = [allRGB; trlCols];                 %#ok<AGROW>
    allDates = [allDates; repmat(string(dateStr), size(trlCols,1), 1)]; %#ok<AGROW>
end

if isempty(allRGB)
    error('Did not collect any DiscprobeColor RGBs from the scanned sessions.');
end

% unique RGBs across all sessions
[uniqRGB, ~, ic] = unique(allRGB, 'rows');
nColors = size(uniqRGB,1);

% how many trials used each RGB
trialCount = accumarray(ic, 1, [nColors 1]);

% earliest date each RGB was seen
firstDate = strings(nColors,1);
for k = 1:nColors
    dks = allDates(ic == k);
    if isempty(dks)
        firstDate(k) = "";
    else
        firstDate(k) = min(dks);
    end
end

% which of these RGBs are already in the legacy ProbeCols table
[isLegacy, legacyRow] = ismember(uniqRGB, legacyRGB, 'rows');

extraRGB      = uniqRGB(~isLegacy,:);
extraCount    = trialCount(~isLegacy);
extraFirst    = firstDate(~isLegacy);

fprintf('\nFound %d unique RGBs total.\n', nColors);
fprintf('  %d already in DiscProbeColsUpdated.\n', sum(isLegacy));
fprintf('  %d new RGBs not in DiscProbeColsUpdated.\n', numel(extraCount));

% pack output struct
META = struct();
META.legacy.ProbeCols   = legacyRGB;
META.legacy.ProbeColIDs = legacyIDs;

META.uniqueRGB      = uniqRGB;
META.uniqueCount    = trialCount;
META.uniqueFirstDate = firstDate;
META.isLegacy       = isLegacy;
META.legacyRow      = legacyRow;

META.extraRGB       = extraRGB;
META.extraCount     = extraCount;
META.extraFirstDate = extraFirst;

% handy table of the new stuff
META.extraTable = table( ...
    extraRGB(:,1), extraRGB(:,2), extraRGB(:,3), ...
    extraCount, extraFirst, ...
    'VariableNames', {'R','G','B','nTrials','firstDate'});

% save to disk for future upgrades
outMat = fullfile(config.paths.code, 'DiscProbeColorMetadata_extended.mat');
outCSV = fullfile(config.paths.code, 'DiscProbe_extra_colors.csv');

try
    save(outMat, 'META');
    fprintf('Saved extended color metadata to %s\n', outMat);
catch ME
    warning('Could not save %s (%s)', outMat, ME.message);
end

try
    writetable(META.extraTable, outCSV);
    fprintf('Saved extra RGB table to %s\n', outCSV);
catch ME
    warning('Could not write %s (%s)', outCSV, ME.message);
end

end
