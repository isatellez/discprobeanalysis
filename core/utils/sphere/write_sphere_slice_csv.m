function write_sphere_slice_csv(R, sessionID, config)
% Write per-slice stats from sphere_slicer into a session CSV.

    if nargin < 3 || isempty(config)
        config = load_config();
    end

    if isempty(R)
        warning('write_sphere_slice_csv: R is empty, nothing to write.');
        return;
    end

    P = get_session_paths(config, sessionID);
    tablesDir = P.tables;
    if ~exist(tablesDir,'dir')
        mkdir(tablesDir);
    end

    n = numel(R);
    sessionID_col = strings(n,1);
    dateStr_col   = strings(n,1);
    monkey_col    = strings(n,1);
    unitIdx_col   = nan(n,1);
    unitID_col    = nan(n,1);
    unitType_col  = strings(n,1);
    sat_col       = nan(n,1);
    elev_col      = nan(n,1);
    nTrials_col   = nan(n,1);
    r2_flat_col   = nan(n,1);
    r2_1_col      = nan(n,1);
    r2_12_col     = nan(n,1);
    bestR2_col    = nan(n,1);
    rayR_col      = nan(n,1);
    rayP_col      = nan(n,1);
    rayMu_col     = nan(n,1);

    for i = 1:n
        sessionID_col(i) = string(R(i).sessionID);
        dateStr_col(i)   = string(R(i).dateStr);
        monkey_col(i)    = string(R(i).monkey);
        unitIdx_col(i)   = R(i).unitIdx;
        unitID_col(i)    = R(i).unitID;
        unitType_col(i)  = string(R(i).unitType);
        sat_col(i)       = R(i).satID;
        elev_col(i)      = R(i).elevID;
        nTrials_col(i)   = R(i).nTrials;

        if isfield(R(i),'r2_flat'),   r2_flat_col(i) = R(i).r2_flat;   end
        if isfield(R(i),'r2_1'),      r2_1_col(i)    = R(i).r2_1;      end
        if isfield(R(i),'r2_12'),     r2_12_col(i)   = R(i).r2_12;     end
        if isfield(R(i),'rayleigh_R'),    rayR_col(i)  = R(i).rayleigh_R;     end
        if isfield(R(i),'rayleigh_p'),    rayP_col(i)  = R(i).rayleigh_p;     end
        if isfield(R(i),'rayleigh_mu_deg'), rayMu_col(i) = R(i).rayleigh_mu_deg; end

        % best-model R2: pick between 1st and 1st+2nd harmonic
        r1  = r2_1_col(i);
        r12 = r2_12_col(i);
        if isfinite(r1) || isfinite(r12)
            if ~isfinite(r12) || (isfinite(r1) && r1 >= r12)
                bestR2_col(i) = r1;
            else
                bestR2_col(i) = r12;
            end
        end
    end

    T = table( ...
        sessionID_col, dateStr_col, monkey_col, ...
        unitIdx_col, unitID_col, unitType_col, ...
        sat_col, elev_col, nTrials_col, ...
        r2_flat_col, r2_1_col, r2_12_col, bestR2_col, ...
        rayR_col, rayP_col, rayMu_col, ...
        'VariableNames', { ...
        'sessionID','dateStr','monkey', ...
        'unitIdx','unitID','unitType', ...
        'satID','elevID','nTrials', ...
        'r2_flat','r2_1','r2_12','bestR2', ...
        'rayleigh_R','rayleigh_p','rayleigh_mu_deg'});

    sessTag = char(sessionID);
    outCSV = fullfile(tablesDir, sprintf('%s_sphereSlices.csv', sessTag));
    writetable(T, outCSV);
    fprintf('write_sphere_slice_csv: wrote %d rows to %s\n', height(T), outCSV);
end
