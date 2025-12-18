function run_sphere_slicer(sessionID, config)
% Sphere Slicer: Generates tuning figures for DiscProbe data.
% SPEED OPTIMIZED: Removed redundant plotting loops.
% STABILITY: Streaming Mode (Writes to disk immediately).

    % --- 1. GLOBAL SILENCE START ---
    origVis = get(0, 'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off');
    set(0, 'DefaultFigureWindowStyle', 'normal');
    cleanupObj = onCleanup(@() set(0, 'DefaultFigureVisible', origVis));

    if nargin < 2 || isempty(config)
        config = load_sphere_config();
    end
    
    % Force the config flag for good measure
    if ~isfield(config, 'plot'), config.plot = struct(); end
    config.plot.makePlots = false; 

    if nargin < 1 || isempty(sessionID)
        sessionID = config.sphere.defaultSession;
    end

    COL = load_color_metadata(config);
    S = load_session(sessionID, config);
    
    if ~isfield(S, 'sessionID'), S.sessionID = sessionID; end
    if ~isfield(S, 'dateStr')
        tok = regexp(sessionID, '\d{6}', 'match', 'once');
        if ~isempty(tok), S.dateStr = tok; else, S.dateStr = sessionID; end
    end

    Tidx = load_or_make_trial_index(S, COL, config);
    nUnits = numel(S.spk);

    % Output Paths
    SessPaths = get_session_paths(config, S.sessionID);
    csvDir = fullfile(SessPaths.tables, 'sphere_slices');
    if ~exist(csvDir, 'dir'), mkdir(csvDir); end
    csvFile = fullfile(csvDir, sprintf('%s_sphere_slices.csv', S.sessionID));

    % --- STREAMING SETUP: Start Fresh ---
    if exist(csvFile, 'file')
        fprintf('[sphere_slicer] Deleting old CSV to start fresh stream: %s\n', csvFile);
        delete(csvFile);
    end

    sessFigRoot = fullfile(SessPaths.figs, 'sphere_slices');
    
    fprintf('[sphere_slicer] Session %s: Found %d units. Starting Stream...\n', S.sessionID, nUnits);

    % --- QC SETTINGS ---
    doOutliers = false;
    if isfield(config, 'qc') && isfield(config.qc, 'removeOutliers')
        doOutliers = logical(config.qc.removeOutliers);
    end
    minTrialsSlice = 10;
    if isfield(config, 'trials') && isfield(config.trials, 'minTrialsSlice')
        minTrialsSlice = config.trials.minTrialsSlice;
    end
    
    % This helper allows 'dklSuite' to default to TRUE if not found in config
    checkFlag = @(field) (~isfield(config, 'sphere') || ~isfield(config.sphere, 'do') || ...
                          ~isfield(config.sphere.do, field) || config.sphere.do.(field));

    % --- UNIT LOOP ---
    for uu = 1:nUnits
        
        % Safe ID printing
        if isfield(S, 'phyIDs'), currID = S.phyIDs(uu); 
        elseif isfield(S, 'unitIDs'), currID = S.unitIDs(uu); 
        else, currID = uu; end

        fprintf('  Processing Unit %d / %d (ID: %d)... ', uu, nUnits, currID);
        
        % Local batch for this unit ONLY
        unit_batch = {}; 
        
        try
            U = prepare_unit(S, uu, Tidx, config);
        catch
            fprintf('Skipped (Load Error)\n');
            continue;
        end
        if doOutliers, U = detect_outliers(U, config); end

        if ~isfield(U, 'trials') || isempty(U.trials) || ...
           ~all(ismember({'hueID','satID','elevID'}, U.trials.Properties.VariableNames))
            fprintf('Skipped (No Trials)\n');
            continue;
        end

        hueID = U.trials.hueID; satID = U.trials.satID; elevID = U.trials.elevID;
        satVals = unique(satID(~isnan(satID))); 
        elevVals = unique(elevID(~isnan(elevID)));

        if isempty(satVals) || isempty(elevVals)
            fprintf('Skipped (No Valid Conditions)\n');
            continue; 
        end
        
        fprintf('OK\n'); 

        unitFolderName = sprintf('%s%d', U.unitType, U.unitID);

        % --- ELEVATION LOOP (PLOTTING) ---
        for ei = 1:numel(elevVals)
            eVal = elevVals(ei);
            usePlane = elevID == eVal & ~isnan(hueID) & ~isnan(satID);
            if nnz(usePlane) < minTrialsSlice, continue; end

            Uplane = U;
            Uplane.spk = U.spk(usePlane); Uplane.mrk = U.mrk(usePlane);
            Uplane.trials = U.trials(usePlane,:); Uplane.nTrials = numel(Uplane.spk);

            % --- SINGLE SAVE DESTINATION ---
            planeSubDir = fullfile('isoluminant_planes', sprintf('elev_%0.2f', eVal));
            saveDir     = fullfile(sessFigRoot, planeSubDir);
            
            if ~exist(saveDir, 'dir'), mkdir(saveDir); end
            
            % Ensure sibling "hues" folder exists (for spike accounting)
            [isoPath, ~] = fileparts(saveDir); 
            [rootPath, ~] = fileparts(isoPath); 
            if ~exist(fullfile(rootPath, 'hues'), 'dir'), mkdir(fullfile(rootPath, 'hues')); end

            plotConfig = config;
            if ~isfield(plotConfig, 'plot'), plotConfig.plot = struct(); end
            plotConfig.plot.makePlots = false; 

            % --- PLOT 1: RAYLEIGH ---
            if checkFlag('rayleighPlot')
                set(0, 'DefaultFigureVisible', 'off');
                try
                    plot_rayleigh(Uplane, struct(), COL, plotConfig, saveDir);
                    close all; 
                catch ME
                    fprintf('  !!! ERROR in plot_rayleigh for Unit %d: %s\n', currID, ME.message);
                    close all;
                end
            end

            % --- PLOT 2: SPIKE ACCOUNTING ---
            if checkFlag('spikeAccounting')
                set(0, 'DefaultFigureVisible', 'off');
                try
                    plot_spike_accounting(Uplane, struct(), COL, plotConfig, saveDir, '');
                    close all; 
                catch ME
                     fprintf('  !!! ERROR in plot_spike_accounting for Unit %d: %s\n', currID, ME.message);
                     close all;
                end
            end

            % --- PLOT 3: DKL SUITE (NEW) ---
            if checkFlag('dklSuite')
                set(0, 'DefaultFigureVisible', 'off');
                try
                    % Plots DKL bars for this specific elevation plane
                    plot_dkl_suite(Uplane, struct(), COL, plotConfig, saveDir);
                    close all; 
                catch ME
                     fprintf('  !!! ERROR in plot_dkl_suite for Unit %d: %s\n', currID, ME.message);
                     close all;
                end
            end
            
            clear Uplane;
        end

        % --- CALC STATS (SLICES) ---
        for si = 1:numel(satVals)
            sVal = satVals(si);
            for ei = 1:numel(elevVals)
                eVal = elevVals(ei);
                useSlice = satID == sVal & elevID == eVal & ~isnan(hueID);
                if nnz(useSlice) < minTrialsSlice, continue; end

                Usub = U;
                Usub.spk = U.spk(useSlice); Usub.mrk = U.mrk(useSlice);
                Usub.trials = U.trials(useSlice,:); Usub.nTrials = numel(Usub.spk);
                
                hm = []; P = []; COS = []; RY = [];

                if checkFlag('hueMeans')
                    try hm = compute_hue_means(Usub); catch, hm = []; end
                end
                
                if isempty(hm) && checkFlag('hueMeans'), continue; end
                
                if checkFlag('peakModel')
                    try P = peak_model_harmonic(Usub, config); catch, P = []; end
                end
                if checkFlag('cosineFit')
                    try COS = cosine_fit(Usub, config); catch, COS = []; end
                end
                if checkFlag('rayleighStats')
                    try RY = rayleigh_test(Usub, config, COL); catch, RY = []; end
                end

                % --- BUILD ROW ---
                row = struct();
                row.sessionID = S.sessionID; row.dateStr = S.dateStr; row.monkey = S.monkey;
                row.unitIdx = U.idx; row.unitID = U.unitID; row.unitType = U.unitType;
                row.satID = sVal; row.elevID = eVal; row.nTrials = Usub.nTrials;
                
                if ~isempty(hm)
                    row.hues = hm.hues; 
                    row.rate_mean = hm.rate_mean;
                    if isfield(hm, 'rate_sem'), row.rate_sem = hm.rate_sem; else, row.rate_sem = []; end
                else
                    row.hues = []; row.rate_mean = []; row.rate_sem = [];
                end

                if ~isempty(P)
                    row.peakClass = P.class;
                    if isfield(P,'p_tuned'), row.p_tuned = P.p_tuned; else, row.p_tuned = NaN; end
                    if isfield(P,'p_bimodal'), row.p_bimodal = P.p_bimodal; else, row.p_bimodal = NaN; end
                    if isfield(P,'r2_flat'), row.r2_flat = P.r2_flat; else, row.r2_flat = NaN; end
                    if isfield(P,'r2_1'), row.r2_1 = P.r2_1; else, row.r2_1 = NaN; end
                    if isfield(P,'r2_12'), row.r2_12 = P.r2_12; else, row.r2_12 = NaN; end
                else
                    row.peakClass = NaN; row.p_tuned = NaN; row.p_bimodal = NaN;
                    row.r2_flat = NaN; row.r2_1 = NaN; row.r2_12 = NaN;
                end

                if ~isempty(COS)
                    if isfield(COS,'pref_deg'), row.cos_pref_deg = COS.pref_deg; else, row.cos_pref_deg = NaN; end
                    if isfield(COS,'amp'), row.cos_amp = COS.amp; else, row.cos_amp = NaN; end
                    if isfield(COS,'R2'), row.cos_R2 = COS.R2; else, row.cos_R2 = NaN; end
                else
                    row.cos_pref_deg = NaN; row.cos_amp = NaN; row.cos_R2 = NaN;
                end

                if ~isempty(RY)
                    if isfield(RY,'R'), row.rayleigh_R = RY.R; else, row.rayleigh_R = NaN; end
                    if isfield(RY,'p'), row.rayleigh_p = RY.p; else, row.rayleigh_p = NaN; end
                    if isfield(RY,'mu_deg'), row.rayleigh_mu = RY.mu_deg; else, row.rayleigh_mu = NaN; end
                else
                    row.rayleigh_R = NaN; row.rayleigh_p = NaN; row.rayleigh_mu = NaN;
                end
                
                unit_batch{end+1} = row;
                clear Usub hm P COS RY;
            end
        end
        
        % --- STREAM TO DISK ---
        if ~isempty(unit_batch)
            try
                T_batch = struct2table([unit_batch{:}]);
                if ~exist(csvFile, 'file')
                    writetable(T_batch, csvFile); 
                else
                    writetable(T_batch, csvFile, 'WriteMode', 'Append'); 
                end
            catch ME
                fprintf('  !!! ERROR writing CSV for Unit %d: %s\n', currID, ME.message);
            end
        end
        
        clear U unit_batch T_batch; 
        if mod(uu, 5) == 0
            java.lang.System.gc();
        end
    end

    fprintf('[sphere_slicer] Done. Slices saved to %s\n', csvFile);
end