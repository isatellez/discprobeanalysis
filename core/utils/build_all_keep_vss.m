function build_all_keep_vss(DATES)
% Build ALL_keep_vss.csv from per-session VSS peakModel CSVs.
%
% Looks for:
%   <DiscProbe base>/vss/<date>/tables/<date>_vss_peakModel_sat1.csv
% and writes:
%   <DiscProbe base>/ALL_keep_vss.csv
%
% Each row = one unit.
% keep = 1 for unimodal/bimodal units, 0 otherwise.

clc;

rootDir = fileparts(mfilename('fullpath'));
addpath(genpath(rootDir));

config = load_config();

if ~isfield(config,'paths') || ~isfield(config.paths,'base') || isempty(config.paths.base)
    error('build_all_keep_vss:NoBasePath', ...
        'config.paths.base is missing or empty; cannot locate vss folder.');
end

baseRoot = config.paths.base;          % e.g. /Users/.../DiscProbe
vssRoot  = fullfile(baseRoot, 'vss');  % e.g. /Users/.../DiscProbe/vss

% handle DATES input
if nargin < 1 || isempty(DATES)
    if ~isfolder(vssRoot)
        warning('build_all_keep_vss:NoVSSRoot', ...
            'No VSS root folder found at %s', vssRoot);
        DATES = {};
    else
        D = dir(vssRoot);
        isDir  = [D.isdir];
        names  = {D.name};
        names  = names(isDir & ~ismember(names, {'.','..'}));
        isDate = cellfun(@(s) ~isempty(regexp(s, '^\d{6}$', 'once')), names);
        DATES  = sort(names(isDate));
    end
else
    if ischar(DATES) || isstring(DATES)
        DATES = {char(DATES)};
    end
end

allKeepRows = table();  % accumulate across dates

for d = 1:numel(DATES)
    dateStr    = char(DATES{d});
    tablesRoot = fullfile(vssRoot, dateStr, 'tables');
    csvPath    = fullfile(tablesRoot, sprintf('%s_vss_peakModel_sat1.csv', dateStr));

    if ~isfile(csvPath)
        warning('build_all_keep_vss:NoCSV', ...
            'No VSS peakModel CSV for %s at %s', dateStr, csvPath);
        continue;
    end

    T = readtable(csvPath);

    if height(T) == 0
        continue;
    end

    % need unitLabel and some kind of peakClass
    if ~ismember('unitLabel', T.Properties.VariableNames)
        warning('build_all_keep_vss:NoUnitLabel', ...
            'Table %s has no unitLabel column, skipping', csvPath);
        continue;
    end

    % prefer the original peak_model classification if present
    if ismember('peakClass_orig', T.Properties.VariableNames)
        cls = string(T.peakClass_orig);
    elseif ismember('peakClass', T.Properties.VariableNames)
        cls = string(T.peakClass);
    else
        warning('build_all_keep_vss:NoPeakClass', ...
            'Table %s has no peakClass/peakClass_orig column, skipping', csvPath);
        continue;
    end

    n = height(T);

    % everyone gets a row; keep = 1 only for uni/bi
    isUniBi = ismember(cls, ["unimodal","bimodal"]);
    keepCol = double(isUniBi);

    dateCol  = repmat(string(dateStr), n, 1);
    unitCol  = string(T.unitLabel(:));

    if ismember('unitType', T.Properties.VariableNames)
        unitTypeCol = string(T.unitType(:));
    else
        unitTypeCol = repmat(missing, n, 1);
    end

    if ismember('unitID', T.Properties.VariableNames)
        unitIDCol = T.unitID(:);
    else
        unitIDCol = nan(n, 1);
    end

    peakClassCol = cls(:);

    keepRowsThisDate = table( ...
        dateCol, ...
        unitCol, ...
        keepCol, ...
        peakClassCol, ...
        unitTypeCol, ...
        unitIDCol, ...
        'VariableNames', {'session','unit','keep','peakClass','unitType','unitID'});

    allKeepRows = [allKeepRows; keepRowsThisDate]; %#ok<AGROW>
end

keepFile = fullfile(baseRoot, 'ALL_keep_vss.csv');

if isempty(allKeepRows)
    allKeepRows = table( ...
        strings(0,1), ...   % session
        strings(0,1), ...   % unit
        zeros(0,1),   ...   % keep
        strings(0,1), ...   % peakClass
        strings(0,1), ...   % unitType
        zeros(0,1),   ...   % unitID
        'VariableNames', {'session','unit','keep','peakClass','unitType','unitID'});
end

try
    writetable(allKeepRows, keepFile);
    fprintf('Wrote ALL_keep_vss.csv with %d rows to:\n  %s\n', ...
        height(allKeepRows), keepFile);
catch ME
    warning('build_all_keep_vss:WriteFailed', ...
        'Could not write ALL_keep_vss.csv: %s', ME.message);
end

end
