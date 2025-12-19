function ColorTable = get_color_coords_table(config)
    % 1. Load IDs
    colsFile = fullfile(config.paths.code, 'DiscProbeColsUpdated.mat');
    load(colsFile, 'ProbeColIDs'); % [Hue, Sat, Elev]
    
    % 2. Load Coordinates
    lmsFile = fullfile(config.paths.code, 'lms.csv');
    LMS_all = readmatrix(lmsFile); % [L, M, S]
    
    % 3. Join them (Assuming rows match 1-to-1)
    % If they don't match, we'd use RGB matching, but usually they are aligned.
    ColorTable = table(ProbeColIDs(:,1), ProbeColIDs(:,2), ProbeColIDs(:,3), ...
                       LMS_all(:,1), LMS_all(:,2), LMS_all(:,3), ...
                       'VariableNames', {'hueID','satVal','elevID','L','M','S'});
end