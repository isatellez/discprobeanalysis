function S = load_session(sessionID, config)
% Load one DiscProbe session and pack spikes into S.
% sessionID can be '250513' or 'Jacomo_250513'.

    if isstring(sessionID)
        sessionID = char(sessionID);
    end

    info = parse_session_id(sessionID); %turns session id into components

    % date for file search
    dateStr = char(info.dateStr);

    % monkey from folder
    % but there is also a fall back to config.monkey for scripts that use old
    % analysis
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

    %we're looking for ExptTrialsDisc and ExptInfo
    if ~isfield(raw, 'ExptTrialsDisc') || ~isfield(raw, 'ExptInfo')
        error('Expected ExptTrialsDisc and ExptInfo in %s', expFile);
    end

    ExptTrialsDisc = raw.ExptTrialsDisc;
    ExptInfo       = raw.ExptInfo;

    if ~isfield(ExptInfo, 'nSU') || ~isfield(ExptInfo, 'nMU')
        error('ExptInfo is missing nSU/nMU in %s', expFile);
    end

    %single units, multi units, total units
    nSU  = ExptInfo.nSU;
    nMU  = ExptInfo.nMU;
    nTOT = nSU + nMU;

    %total trials
    nTrials = size(ExptTrialsDisc, 1);
    if nTOT <= 0 || nTrials == 0
        error('No units or no trials in %s', expFile);
    end


    cols = 10:(9+nTOT); %theres a the cols that have spike times
    if size(ExptTrialsDisc, 2) < max(cols)
        error('Expected at least %d columns in ExptTrialsDisc, got %d.', ...
            max(cols), size(ExptTrialsDisc,2));
    end

    if isfield(config, 'time') && isfield(config.time, 'fullWin')
        winFull = config.time.fullWin;
    else
        winFull = [0 0.4]; %this is the real time window
    end


    spkCols = ExptTrialsDisc(:, cols); %spikes
    trimfun = @(x) x(x >= winFull(1) & x <= winFull(2)); %trims to desired window
    spkCols = cellfun(trimfun, spkCols, 'UniformOutput', false);

    spk_all = cell(1, nTOT);
    mrk_all = cell(1, nTOT);

    mrkWin   = [winFull(1) winFull(2)]; 
    mrk_cell = repmat({mrkWin}, nTrials, 1);  %stack of windows for rasters

    for uu = 1:nTOT %for every unit
        unit_tr = spkCols(:, uu);
        for tt = 1:nTrials %for every trial 
            v = unit_tr{tt}; 
            if isempty(v)
                unit_tr{tt} = [];
            else
                unit_tr{tt} = double(v(:).'); %row of spikes 
            end
        end
        spk_all{uu} = unit_tr; %stash spikes in unit
        mrk_all{uu} = mrk_cell;
    end

    if isfield(ExptInfo, 'spk_ID_SU') && isfield(ExptInfo, 'spk_ID_MU')
        unitIDs_all = [ExptInfo.spk_ID_SU(:); ExptInfo.spk_ID_MU(:)];
    else
        unitIDs_all = (1:nTOT).';
    end

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

    S.phyIDs        = unitIDs_all;
    S.trNum        = (1:nTrials).';
end

