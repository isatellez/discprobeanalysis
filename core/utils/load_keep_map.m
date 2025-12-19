function keepMap = load_keep_map(baseDiscProbeLocal, primaryName, fallbackName)
% LOAD_KEEP_MAP - Load keep-list CSV into a containers.Map
% Key: 'Session_UnitLabel' (e.g., 'Jacomo_250513_SU1' or '250513_SU1')
% Val: 1 (Keep) or 0 (Drop)
%
% Auto-searches baseDiscProbeLocal AND ../output/tables if not found.

if nargin < 2 || isempty(primaryName)
    primaryName = 'ALL_keep_vss.csv';
end
if nargin < 3 || isempty(fallbackName)
    fallbackName = 'ALL_keep.csv';
end

keepMap = containers.Map;

% Define Candidate Paths
% 1. Direct base path (data/)
p1 = fullfile(baseDiscProbeLocal, primaryName);
p2 = fullfile(baseDiscProbeLocal, fallbackName);

% 2. Output/Tables path (assuming standard ../output/tables structure relative to code or base)
% Try to deduce output/tables from the base path hierarchy or common config
% (Simple heuristic: look in sibling 'output/tables' or subfolder 'tables')
p3 = fullfile(baseDiscProbeLocal, '..', 'output', 'tables', primaryName);
p4 = fullfile(baseDiscProbeLocal, 'tables', primaryName);

if isfile(p1)
    keepFile = p1;
elseif isfile(p3) % Check output folder
    keepFile = p3;
elseif isfile(p4)
    keepFile = p4;
elseif isfile(p2) % Fallback name
    keepFile = p2;
    warning('load_keep_map:Fallback', 'Primary keep list not found. Using fallback: %s', p2);
else
    warning('load_keep_map:NotFound', ...
        'Could not find keep list (%s) in %s or standard output folders. Processing ALL units.', ...
        primaryName, baseDiscProbeLocal);
    return;
end

fprintf('Loading keep list from: %s\n', keepFile);
T_keep = readtable(keepFile);
names      = T_keep.Properties.VariableNames;
namesLower = lower(names);

% --- Detect Columns ---
sessField = '';
for cand = {'session','date','datestr'}
    idx = find(strcmp(namesLower, cand{1}), 1);
    if ~isempty(idx), sessField = names{idx}; break; end
end
if isempty(sessField), error('Keep list missing session/date column.'); end

unitField = '';
for cand = {'unit','unit_label','unitname','unitlabel'}
    idx = find(strcmp(namesLower, cand{1}), 1);
    if ~isempty(idx), unitField = names{idx}; break; end
end
if isempty(unitField), error('Keep list missing unit column.'); end

keepField = '';
for cand = {'keep_vss','keep','keep_flag'}
    idx = find(strcmp(namesLower, cand{1}), 1);
    if ~isempty(idx), keepField = names{idx}; break; end
end
if isempty(keepField), error('Keep list missing keep flag column.'); end

% --- Populate Map ---
for i = 1:height(T_keep)
    % 1. Parse Unit Name
    uName = string(T_keep.(unitField)(i));
    
    % 2. Parse Session String (Handle 250513 vs 'Jacomo_250513')
    sVal = T_keep.(sessField)(i);
    if isnumeric(sVal) && ~isnan(sVal)
        sessStr = sprintf('%06d', sVal);
    else
        sessStr = char(string(sVal));
    end
    
    % 3. Parse Keep Value (Ensure Numeric 1/0)
    rawKeep = T_keep.(keepField)(i);
    if iscell(rawKeep) || isstring(rawKeep) || ischar(rawKeep)
        % Handle "TRUE"/"FALSE" strings
        strK = string(rawKeep);
        if strcmpi(strK, "true") || strcmpi(strK, "1")
            val = 1;
        else
            val = 0;
        end
    else
        val = double(rawKeep);
    end
    
    % Generate Key
    key = sprintf('%s_%s', sessStr, uName);
    keepMap(key) = val;
end
end