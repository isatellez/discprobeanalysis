function ColorTable = get_standardized_color_table(config)
% Master mapping of IDs to physical coordinates.
% Used for joining with trial data or creating 3D model visualizations.

    colsFile = fullfile(config.paths.code, 'DiscProbeColsUpdated.mat');
    load(colsFile, 'ProbeColIDs'); % [Hue, Sat, Elev]
    
    lmsFile = fullfile(config.paths.code, 'lms.csv');
    LMS_all = readmatrix(lmsFile); % [L, M, S]
    
    % Force Standard Names: hueID, satVal, elevID
    ColorTable = table(ProbeColIDs(:,1), ProbeColIDs(:,2), ProbeColIDs(:,3), ...
                       LMS_all(:,1), LMS_all(:,2), LMS_all(:,3), ...
                       'VariableNames', {'hueID','satVal','elevID','L','M','S'});
end