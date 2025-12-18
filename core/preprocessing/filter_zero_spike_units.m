function S = filter_zero_spike_units(S, config)
% drop units with 0 spikes across all trials in the full window

if ~isfield(S, 'spk') || isempty(S.spk)
    return;
end

if isfield(config, 'time') && isfield(config.time, 'fullWin')
    winFull = config.time.fullWin;
else
    winFull = [0 0.4]; %default time window
end

t0  = winFull(1);
t1  = winFull(2);
nU  = S.nUnits;
nTr = S.nTrials;

totalPerUnit = zeros(1, nU);

for uu = 1:nU %for every unit
    spk_cell = S.spk{uu};
    if numel(spk_cell) ~= nTr %does the unit have n trials?
        error('S.spk{%d} has %d trials, expected %d.', uu, numel(spk_cell), nTr);
    end

    c = 0;
    for tt = 1:nTr %for every trial
        spks = spk_cell{tt};
        if isempty(spks)
            continue;
        end
        spks = spks(spks >= t0 & spks < t1);
        c = c + numel(spks);
    end
    totalPerUnit(uu) = c;
end

keepUnits = totalPerUnit > 0; %we keep all units whose total is greater than 0

if all(keepUnits)
    return;
end

fprintf('Removing %d/%d units with 0 spikes.\n', sum(~keepUnits), nU);

S.spk= S.spk(keepUnits);
if isfield(S, 'mrk') && numel(S.mrk) == nU
    S.mrk = S.mrk(keepUnits);
end
if isfield(S, 'unitIDs') && numel(S.phyIDs) == nU
    S.phyIDs = S.phyIDs(keepUnits);
end

S.nUnits = sum(keepUnits);

if isfield(S, 'ExptInfo')
    EI = S.ExptInfo;

    %aligning nsu and msu so everyone is in the same page
    if isfield(EI, 'nSU') && isfield(EI, 'nMU') && ...
       isfield(EI, 'spk_ID_SU') && isfield(EI, 'spk_ID_MU')

        nSU = EI.nSU;
        nMU = EI.nMU;

        if numel(keepUnits) == nSU + nMU
            suMask = keepUnits(1:nSU);
            muMask = keepUnits(nSU+1:nSU+nMU);

            EI.spk_ID_SU = EI.spk_ID_SU(suMask);
            EI.spk_ID_MU = EI.spk_ID_MU(muMask);

            EI.nSU = numel(EI.spk_ID_SU);
            EI.nMU = numel(EI.spk_ID_MU);
        end
    end

    EI.unit_new2old = find(keepUnits);
    S.ExptInfo      = EI;
end

if S.nUnits == 0
    error('All units had zero spikes after filtering. Nothing to analyze.');
end

end
