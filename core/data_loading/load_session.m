function S = load_session(sessionID, config)
% LOAD_SESSION - Load one DiscProbe session and pack spikes into S.
% S.phyIDs    -> The specific Unit IDs (e.g. 13, 42)
% S.unitTypes -> Cell array {'SU', 'MU'} matching the IDs

    if isstring(sessionID)
        sessionID = char(sessionID);
    end
    info = parse_session_id(sessionID); 
    
    dateStr = char(info.dateStr);
    
    if strlength(info.monkey) > 0
        monkey = char(info.monkey);
    else
        monkey = config.monkey;
    end
    
    PATHS = struct();
    PATHS.baseDiscProbeLocal = config.paths.base;
    PATHS.baseDiscProbesCode = config.paths.code;
    PATHS.baseOut            = config.paths.output;
    
    expFile = locate_expTrialsDisc(PATHS, monkey, dateStr);
    if ~isfile(expFile)
        error('Could not find ExpTrialsDisc file for %s (got: %s)', sessionID, expFile);
    end
    
    t = tic;
    raw = load(expFile);
    fprintf('Loaded %s in %.1f s.\n', expFile, toc(t));
    
    if ~isfield(raw, 'ExptTrialsDisc') || ~isfield(raw, 'ExptInfo')
        error('Expected ExptTrialsDisc and ExptInfo in %s', expFile);
    end
    
    ExptTrialsDisc = raw.ExptTrialsDisc;
    ExptInfo       = raw.ExptInfo;
    
    if ~isfield(ExptInfo, 'nSU') || ~isfield(ExptInfo, 'nMU')
        error('ExptInfo is missing nSU/nMU in %s', expFile);
    end
    
    nSU  = ExptInfo.nSU;
    nMU  = ExptInfo.nMU;
    nTOT = nSU + nMU;
    
    nTrials = size(ExptTrialsDisc, 1);
    if nTOT <= 0 || nTrials == 0
        error('No units or no trials in %s', expFile);
    end
    
    cols = 10:(9+nTOT); 
    if size(ExptTrialsDisc, 2) < max(cols)
        error('Expected at least %d columns in ExptTrialsDisc, got %d.', ...
            max(cols), size(ExptTrialsDisc,2));
    end
    
    if isfield(config, 'time') && isfield(config.time, 'fullWin')
        winFull = config.time.fullWin;
    else
        winFull = [0 0.4]; 
    end
    
    spkCols = ExptTrialsDisc(:, cols); 
    trimfun = @(x) x(x >= winFull(1) & x <= winFull(2)); 
    spkCols = cellfun(trimfun, spkCols, 'UniformOutput', false);
    
    spk_all = cell(1, nTOT);
    mrk_all = cell(1, nTOT);
    mrkWin   = [winFull(1) winFull(2)]; 
    mrk_cell = repmat({mrkWin}, nTrials, 1); 
    
    for uu = 1:nTOT 
        unit_tr = spkCols(:, uu);
        for tt = 1:nTrials 
            v = unit_tr{tt}; 
            if isempty(v)
                unit_tr{tt} = [];
            else
                unit_tr{tt} = double(v(:).'); 
            end
        end
        spk_all{uu} = unit_tr; 
        mrk_all{uu} = mrk_cell;
    end
    
    % --- CONSTRUCT IDs AND TYPES ---
    if isfield(ExptInfo, 'spk_ID_SU') && isfield(ExptInfo, 'spk_ID_MU')
        phyIDs_all = [ExptInfo.spk_ID_SU(:); ExptInfo.spk_ID_MU(:)];
    else
        phyIDs_all = (1:nTOT).';
    end
    
    % Create unitTypes manually since it doesn't exist in raw file
    unitTypes_all = [repmat({'SU'}, nSU, 1); repmat({'MU'}, nMU, 1)];

    % Output Structure
    S = struct();
    S.sessionID      = info.sessionID;  
    S.dateStr        = dateStr;
    S.monkey         = monkey;
    S.expFile        = expFile;
    S.ExptInfo       = ExptInfo;
    S.ExptTrialsDisc = ExptTrialsDisc;
    S.nTrials        = nTrials;
    S.nUnits         = nTOT;
    S.spk            = spk_all;
    S.mrk            = mrk_all;
    
    % Use phyIDs as the main identifier
    S.phyIDs         = phyIDs_all;    
    
    % Keep unitIDs as a copy for compatibility
    S.unitIDs        = phyIDs_all;    
    
    % Crucial: Add unitTypes so downstream scripts don't fail
    S.unitTypes      = unitTypes_all;  
    
    S.trNum          = (1:nTrials).';
end