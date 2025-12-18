function U = attach_color_index(U, COL)
% Attach U.colIdx = row in COL.probeIDs / COL.lms per trial.
% Assumes U.trials has hueID, saturID, elevID that match COL.probeIDs columns.

nT = height(U.trials);
colIdx = nan(nT,1);

hue  = U.trials.hueID;
sat  = U.trials.satID;
elev = [];

if isfield(U.trials,'elevID')
    elev = U.trials.elevID;
end

ids = COL.probeIDs;  % [144 x ?], assume columns: 1=hue, 2=sat, 5=elev

for t = 1:nT
    if isnan(hue(t)) || isnan(sat(t))
        continue
    end

    match = ids(:,1) == hue(t) & ids(:,2) == sat(t);

    if ~isempty(elev)
        match = match & ids(:,5) == elev(t);
    end

    k = find(match,1,'first');
    if ~isempty(k)
        colIdx(t) = k;
    end
end

U.colIdx = colIdx;
end
