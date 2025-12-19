function DiscProbe_3D_Coordinate_Regression(DATES)
% Fits a 3D Poisson GLM to all trials using LMS coordinates.
% This identifies the 3D sensitivity vector (L, M, S) using the full 
% response manifold across all elevations and saturations.
%
% Standardized version: 
%   Uses core/utils/normalize_dates and resolve_session_folder
%   Uses standardized make_trial_index_table (L, M, S columns)
%   Employs Poisson regression for robust neural counting
    clc;
    rootDir = fileparts(mfilename('fullpath'));
    addpath(genpath(rootDir));
    
    config = load_config();
    COL    = load_color_metadata(config);
    
    % Use standardized date handling
    DATES = normalize_dates(DATES, config); 
    % Analysis window (0.05 to 0.20s typical for stimulus-driven response)
    timeWin = [0.05 0.20]; 
    
    allUnitsResults = [];
    for d = 1:numel(DATES)
        % Use standardized folder resolution
        sessionFolder = resolve_session_folder(config.paths.base, DATES{d});
        
        if isempty(sessionFolder)
            warning('Could not find folder for %s. Skipping.', DATES{d});
            continue;
        end
        
        fprintf('\n=== 3D Trial-by-Trial Analysis: %s ===\n', sessionFolder);
        % Load Session Data
        try
            S = load_session(sessionFolder, config);
        catch ME
            fprintf('  Load failed: %s\n', ME.message);
            continue;
        end
        
        % Generate the Standardized Trial Table
        % This table has T.L, T.M, T.S and T.elevID
        T = make_trial_index_table(S, COL, config);
        
        for cc = 1:S.nUnits
            phyID = S.unitIDs(cc);
            uType = S.unitTypes{cc};
            unitLabel = sprintf('%s%d', uType, phyID);
            
            % Count spikes per trial in the window
            y = cellfun(@(x) sum(x >= timeWin(1) & x <= timeWin(2)), S.spk{cc});
            y = double(y(:)); 
            
            % Predictor Matrix: Standardized LMS coordinates from T
            X = [T.L, T.M, T.S];
            validRows = all(isfinite(X), 2);
            
            % Identify isoluminant subset for internal comparison
            isIso = abs(T.elevID) < 0.1;
            
            X_all = X(validRows, :);
            y_all = y(validRows);
            X_iso = X(validRows & isIso, :);
            y_iso = y(validRows & isIso);

            % filter out units that are too sparse for regression
            if numel(y_all) < 50 || mean(y_all) < 0.1
                continue;
            end

            try
                % Fit using Poisson Identity Link
                % This models the linear summation of cones directly to spike counts
                opts = statset('glmfit');
                opts.MaxIter = 200;
                
                % Full 3D Trial Model
                [b_all, dev_all] = glmfit(X_all, y_all, 'poisson', 'link', 'identity', 'options', opts);
                [~, null_all]    = glmfit(ones(size(y_all,1),1), y_all, 'poisson', 'link', 'identity', 'constant', 'off');
                R2_3D = 1 - (dev_all / null_all);
                
                % Isoluminant Subset Trial Model
                [~, dev_iso]  = glmfit(X_iso, y_iso, 'poisson', 'link', 'identity', 'options', opts);
                [~, null_iso] = glmfit(ones(size(y_iso,1),1), y_iso, 'poisson', 'link', 'identity', 'constant', 'off');
                R2_Iso = 1 - (dev_iso / null_iso);
                
                % report results: R2 values will be lower than mean-fitting but more honest
                fprintf('  %s | R2-3D: %.3f | R2-Iso: %.3f | L:%.2f M:%.2f S:%.2f\n', ...
                    unitLabel, R2_3D, R2_Iso, b_all(2), b_all(3), b_all(4));
                
                allUnitsResults = [allUnitsResults; struct(...
                    'session', sessionFolder, 'unit', unitLabel, ...
                    'R2_3D', R2_3D, 'R2_Iso', R2_Iso, 'bL', b_all(2), 'bM', b_all(3), 'bS', b_all(4))]; %#ok<AGROW>
            catch
                % print failure message in the requested format
                fprintf('  %s | GLM fit failed.\n', unitLabel);
            end
        end
    end
    
    % Save Global Results
    if ~isempty(allUnitsResults)
        outDir = fullfile(config.paths.output, 'tables');
        if ~exist(outDir, 'dir'), mkdir(outDir); end
        outPath = fullfile(outDir, 'ALL_3D_Trial_Regression.csv');
        writetable(struct2table(allUnitsResults), outPath);
        fprintf('\n✅ Finished. Results saved to %s\n', outPath);
    end
end