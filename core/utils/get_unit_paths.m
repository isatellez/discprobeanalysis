function Upaths = get_unit_paths(config, dateStr, unitType, unitID)
    P = get_session_paths(config, dateStr);
    label        = sprintf('%s%d', unitType, unitID);  % "SU10"
    Upaths.root  = fullfile(P.units, label);
    Upaths.mats  = fullfile(Upaths.root, 'mats');
    Upaths.figures = fullfile(Upaths.root, 'figures');
    Upaths.tables  = fullfile(Upaths.root, 'tables');
end
