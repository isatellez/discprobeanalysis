function P = get_session_paths(config, sessionID)
% Build canonical paths for a DiscProbe session.
%
% sessionID can be:
%   '250513'
%   'Jacomo_250513'
%   'Jacomo_250513_am'
%
% UPDATE: Auto-detects 'Monkey_Date' folder if only 'Date' is provided.

    info = parse_session_id(sessionID);
    base = config.paths.base;
    
    sessDir = char(info.sessionID); 

    % --- SMART FOLDER DETECTION ---
    % 1. Check if the folder simply exists as-is (e.g. '250513')
    pathAsIs = fullfile(base, sessDir);
    
    if ~exist(pathAsIs, 'dir')
        % 2. If not, and the ID looks like just a date (no underscore),
        %    search for a matching "Monkey_Date" folder.
        if ~contains(sessDir, '_')
            % Look for any folder ending in this date string (e.g. '*_250513')
            candidates = dir(fullfile(base, ['*_' sessDir]));
            candidates = candidates([candidates.isdir]);
            
            % Remove '.' and '..'
            candidates = candidates(~ismember({candidates.name}, {'.','..'}));
            
            if ~isempty(candidates)
                % Success! Use the existing Monkey_Date folder (e.g. 'Jacomo_250513')
                sessDir = candidates(1).name;
                
                % Update info to match the found folder
                info = parse_session_id(sessDir); 
                
            elseif isfield(config, 'monkey') && ~isempty(config.monkey)
                % 3. Folder doesn't exist, but we have a preferred Monkey in config.
                %    Construct the proper name for the new folder.
                sessDir = sprintf('%s_%s', config.monkey, sessDir);
            end
        end
    end
    % ------------------------------

    P = struct();

    % root session folder
    P.session = fullfile(base, sessDir);

    % standard subfolders
    P.raw     = fullfile(P.session, 'raw');
    P.units   = fullfile(P.session, 'units');
    P.figs    = fullfile(P.session, 'figs');
    P.tables  = fullfile(P.session, 'tables');
    P.stas    = fullfile(P.session, 'stas');
    P.sta_lms = fullfile(P.session, 'sta_lmschannel');

    % expose parsed metadata for convenience
    P.sessionID = sessDir; % Update this to the actual folder name used
    P.dateStr   = info.dateStr;
    P.monkey    = info.monkey;
end