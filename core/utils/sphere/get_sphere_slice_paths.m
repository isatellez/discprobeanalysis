function paths = get_sphere_slice_paths(config, sessionID)
    % session-aware roots
    SPaths = get_session_paths(config, sessionID);

    % 1) under the session root
    paths.sessionRoot = fullfile(SPaths.session, 'sphere_slices');

    % 2) under figs
    paths.figRoot     = fullfile(SPaths.figs, 'sphere_slices');

    % 3) under global output
    if isfield(config, 'paths') && isfield(config.paths, 'output') ...
            && ~isempty(config.paths.output)
        paths.globalRoot = fullfile(config.paths.output, 'sphere_slices', char(SPaths.sessionID));
    else
        paths.globalRoot = '';
    end

    roots = {paths.sessionRoot, paths.figRoot, paths.globalRoot};
    for i = 1:numel(roots)
        if ~isempty(roots{i}) && ~exist(roots{i}, 'dir')
            mkdir(roots{i});
        end
    end
end
