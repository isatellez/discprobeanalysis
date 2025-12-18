function save_vss_unit(config, U, hm, P)
% Save per-unit VSS mat into:
%   <base>/<date>/units/<unitLabel>/mats/vss/VSSUnit_<date>_<unitType>_<unitID>.mat

    dateStr  = char(string(U.dateStr));
    unitType = char(string(U.unitType));
    unitID   = U.unitID;

    % use your helper to get unit paths
    Upaths = get_unit_paths(config, dateStr, unitType, unitID);

    % put VSS mats in a vss subfolder under mats
    vssMatDir = fullfile(Upaths.mats, 'vss');
    if ~exist(vssMatDir, 'dir')
        mkdir(vssMatDir);
    end

    fileTag = sprintf('%s_%s_%d', dateStr, unitType, unitID);
    matFile = fullfile(vssMatDir, sprintf('VSSUnit_%s.mat', fileTag));

    meta = struct();
    meta.dateStr  = dateStr;
    meta.unitType = unitType;
    meta.unitID   = unitID;
    meta.nTrials  = U.nTrials;
    if isfield(U, 'winEarly')
        meta.winEarly = U.winEarly;
    end

    save(matFile, 'U', 'hm', 'P', 'meta');
end
