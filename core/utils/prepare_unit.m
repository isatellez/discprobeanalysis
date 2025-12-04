function U = prepare_unit(S, unitIdx, Tidx, config)
% Build per-unit struct U from session S and trial table Tidx.

U = struct();

% basic metadata
U.dateStr = S.dateStr;
U.unitID  = S.unitIDs(unitIdx);
U.idx     = unitIdx;

% SU vs MU based on ExptInfo counts
if isfield(S, 'ExptInfo') && isfield(S.ExptInfo, 'nSU')
    if unitIdx <= S.ExptInfo.nSU
        U.unitType = 'SU';
    else
        U.unitType = 'MU';
    end
else
    U.unitType = 'unit';
end

% spikes and markers for this unit
spk_full = S.spk{unitIdx};   % cell per trial
mrk_full = S.mrk{unitIdx};   % same shape

nTrials_full = numel(spk_full);

if height(Tidx) ~= nTrials_full
    error('prepare_unit: trial count mismatch (Tidx has %d rows, spikes have %d trials).', ...
        height(Tidx), nTrials_full);
end

% default: keep only trials with valid color IDs
% (matches keep = ~isnan(trType_full(:,1)) in the original v1figproc script)
if ismember('hueID', Tidx.Properties.VariableNames)
    keep = ~isnan(Tidx.hueID);
else
    keep = true(nTrials_full,1);
end

if ~any(keep)
    error('prepare_unit: no trials left after applying color-ID mask.');
end

spk_cell = spk_full(keep);
mrk_cell = mrk_full(keep);
Tunit    = Tidx(keep,:);

% fill or overwrite unitID column for this unit
if ismember('unitID', Tunit.Properties.VariableNames)
    Tunit.unitID(:) = U.unitID;
else
    Tunit.unitID = repmat(U.unitID, height(Tunit), 1);
end

U.nTrials = numel(spk_cell);
U.spk     = spk_cell;
U.mrk     = mrk_cell;
U.trials  = Tunit;

% time windows: pull from config if present, otherwise fall back to sane defaults
if isfield(config, 'time')
    if isfield(config.time, 'winEarly')
        U.winEarly = config.time.winEarly;
    else
        U.winEarly = [0 0.2];
    end

    if isfield(config.time, 'winLate')
        U.winLate = config.time.winLate;
    else
        U.winLate = [0.25 0.4];
    end

    if isfield(config.time, 'fullWin')
        U.winFull = config.time.fullWin;
    else
        U.winFull = [0 0.4];
    end
else
    U.winEarly = [0 0.2];
    U.winLate  = [0.25 0.4];
    U.winFull  = [0 0.4];
end

end
