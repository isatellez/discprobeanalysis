function keepMap = load_keep_map(baseDiscProbeLocal, primaryName, fallbackName)
% Load a keep-list CSV and return a containers.Map keyed by 'YYMMDD_unitLabel'.
%
% baseDiscProbeLocal: root DiscProbe folder (config.paths.base).
% primaryName:  optional, default 'ALL_keep_vss.csv'.
% fallbackName: optional, default 'ALL_keep.csv'.

if nargin < 2 || isempty(primaryName)
    primaryName = 'ALL_keep_vss.csv';
end
if nargin < 3 || isempty(fallbackName)
    fallbackName = 'ALL_keep.csv';
end

keepMap = containers.Map;

keepFilePrimary  = fullfile(baseDiscProbeLocal, primaryName);
keepFileFallback = fullfile(baseDiscProbeLocal, fallbackName);

if isfile(keepFilePrimary)
    keepFile = keepFilePrimary;
elseif isfile(keepFileFallback)
    keepFile = keepFileFallback;
    warning('load_keep_map:Fallback', ...
        'Primary keep list %s not found, using fallback %s instead.', ...
        keepFilePrimary, keepFileFallback);
else
    warning('load_keep_map:NotFound', ...
        'No keep list found (%s or %s). Processing all units.', ...
        keepFilePrimary, keepFileFallback);
    return;
end

fprintf('Loading keep list from %s...\n', keepFile);

T_keep = readtable(keepFile);

names      = T_keep.Properties.VariableNames;
namesLower = lower(names);

% session/date column
sessField = '';
for cand = {'session','date','datestr'}
    idx = find(strcmp(namesLower, cand{1}), 1);
    if ~isempty(idx)
        sessField = names{idx};
        break;
    end
end
if isempty(sessField)
    error('Keep list %s is missing a session/date column.', keepFile);
end

% unit label column
unitField = '';
for cand = {'unit','unit_label','unitname'}
    idx = find(strcmp(namesLower, cand{1}), 1);
    if ~isempty(idx)
        unitField = names{idx};
        break;
    end
end
if isempty(unitField)
    error('Keep list %s is missing a unit column.', keepFile);
end

% keep flag column (prefer VSS-specific if present)
keepField = '';
for cand = {'keep_vss','keep','keep_flag'}
    idx = find(strcmp(namesLower, cand{1}), 1);
    if ~isempty(idx)
        keepField = names{idx};
        break;
    end
end
if isempty(keepField)
    error('Keep list %s is missing a keep flag column.', keepFile);
end

for i = 1:height(T_keep)
    uName = string(T_keep.(unitField)(i));

    sVal = T_keep.(sessField)(i);
    if isnumeric(sVal)
        sessStr = sprintf('%06d', sVal);
    else
        sessStr = char(string(sVal));
    end

    key = sprintf('%s_%s', sessStr, uName);
    keepMap(key) = T_keep.(keepField)(i);
end
end
