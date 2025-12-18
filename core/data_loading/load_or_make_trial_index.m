function Tidx = load_or_make_trial_index(S, COL, config)
% Load or build the per-session TrialIndex table.
%
% If a CSV already exists for this session, load it.
% Otherwise, call make_trial_index_table and save a new CSV.
%

if nargin < 3
    config = struct();
end

% figure out the session name (Jocamo_250513 style if present)
if isfield(S, 'sessionID') && ~isempty(S.sessionID)
    sessName = char(S.sessionID);
else
    sessName = char(S.dateStr);
end

% ---------- figure out where the CSV should live ----------

idxDir = '';

% preferred: use per-session root dir via get_session_paths
if exist('get_session_paths','file') == 2
    try
        P = get_session_paths(config, sessName);
        if isfield(P, 'session') && ~isempty(P.session)
            idxDir = P.session;   % <DiscProbe base>/Jocamo_250513
        end
    catch ME
        warning('load_or_make_trial_index: get_session_paths failed (%s). Falling back.', ME.message);
    end
end

% fall back to base/<session> or output/<session> if needed
if isempty(idxDir)
    if isfield(config, 'paths') && isfield(config.paths, 'base') ...
            && ~isempty(config.paths.base)
        idxDir = fullfile(config.paths.base, sessName);
    elseif isfield(config, 'paths') && isfield(config.paths, 'output') ...
            && ~isempty(config.paths.output)
        idxDir = fullfile(config.paths.output, sessName);
    else
        % no place to save → build in memory only
        Tidx = make_trial_index_table(S, COL, config);
        return;
    end
end

if ~exist(idxDir, 'dir')
    mkdir(idxDir);
end

idxCSV = fullfile(idxDir, sprintf('%s_TrialIndex.csv', sessName));

% ---------- decide whether to reuse or rebuild ----------

forceRebuild = false;
if isfield(config, 'trials') && isfield(config.trials, 'forceRebuildIndex') ...
        && ~isempty(config.trials.forceRebuildIndex)
    forceRebuild = logical(config.trials.forceRebuildIndex);
end

% try to load existing CSV
if isfile(idxCSV) && ~forceRebuild
    try
        Tidx = readtable(idxCSV);
        return;
    catch ME
        warning('load_or_make_trial_index: failed to read %s (%s). Rebuilding.', ...
                idxCSV, ME.message);
    end
end

% build fresh and write it where we just decided
Tidx = make_trial_index_table(S, COL, config);

try
    writetable(Tidx, idxCSV);
catch ME
    warning('load_or_make_trial_index: could not write %s (%s).', idxCSV, ME.message);
end

end
