function save_unit_outputs(config, U, stats)
% Save per-unit DiscProbe outputs into:
%   <base>/<date>/units/<unitLabel>/mats/discprobe/DiscProbeUnit_<date>_<unitType>_<unitID>.mat
% Also optionally append a quick CSV row into tablesRoot (if set in config).

    if nargin < 3
        error('save_unit_outputs expects (config, U, stats).');
    end

    if ~isfield(U, 'unitID') || ~isfield(U, 'unitType') || ~isfield(U, 'dateStr')
        error('U is missing unitID, unitType or dateStr.');
    end

    dateStr  = char(string(U.dateStr));
    unitType = char(string(U.unitType));
    unitID   = U.unitID;

    % per-unit paths via helper
    Upaths = get_unit_paths(config, dateStr, unitType, unitID);

    % put DiscProbe mats under mats/discprobe
    discMatDir = fullfile(Upaths.mats, 'discprobe');
    if ~exist(discMatDir, 'dir')
        mkdir(discMatDir);
    end

    unitLabel = sprintf('%s%d', unitType, unitID);  % e.g. SU_10
    fileTag   = sprintf('%s_%s_%d', dateStr, unitType, unitID);
    matFile   = fullfile(discMatDir, sprintf('DiscProbeUnit_%s.mat', fileTag));

    meta = struct();
    meta.dateStr   = string(U.dateStr);
    meta.unitID    = U.unitID;
    meta.unitType  = string(U.unitType);
    if isfield(U,'nTrials'),  meta.nTrials  = U.nTrials;  end
    if isfield(U,'winEarly'), meta.winEarly = U.winEarly; end
    if isfield(U,'winLate'),  meta.winLate  = U.winLate;  end
    if isfield(U,'winFull'),  meta.winFull  = U.winFull;  end
    if isfield(U,'phyID'),    meta.phyID    = U.phyID;    end

    try
        save(matFile, 'U', 'stats', 'meta', '-v7');
    catch ME
        warning('Could not save unit outputs to %s: %s', matFile, ME.message);
    end

    % ---------- optional quick CSV summary ----------
    % we'll now use config.paths.tablesRoot (set in run_discprobe_analysis)
    tablesRoot = '';
    if isfield(config,'paths') && isfield(config.paths,'tablesRoot') ...
            && ~isempty(config.paths.tablesRoot)
        tablesRoot = config.paths.tablesRoot;
    end

    if isempty(tablesRoot)
        return;
    end

    if ~exist(tablesRoot, 'dir')
        mkdir(tablesRoot);
    end

    try
        csvFile    = fullfile(tablesRoot, sprintf('%s_unit_quickSummary.csv', dateStr));
        writeHeader = ~isfile(csvFile);

        fid = fopen(csvFile, 'a');
        if fid == -1
            error('Could not open %s for writing.', csvFile);
        end

        if writeHeader
            hdr = [
                "dateStr", "unitType", "unitID", "nTrials", ...
                "meanRateOverall", ...
                "R_rayleigh", "p_rayleigh", ...
                "R2_cosine", ...
                "p_perm", ...
                "p_anova_hue", "p_anova_sat", "p_anova_int", ...
                "peakClass" ...
            ];
            fprintf(fid, '%s\n', strjoin(hdr, ","));
        end

        meanRateOverall = NaN;
        if isfield(stats, 'hueMeans') && isfield(stats.hueMeans, 'rate_mean')
            meanRateOverall = mean(stats.hueMeans.rate_mean, 'omitnan');
        end

        R_rayleigh = NaN;
        p_rayleigh = NaN;
        if isfield(stats, 'rayleigh') && ~isempty(stats.rayleigh)
            if isfield(stats.rayleigh, 'R'), R_rayleigh = stats.rayleigh.R; end
            if isfield(stats.rayleigh, 'p'), p_rayleigh = stats.rayleigh.p; end
        end

        R2_cosine = NaN;
        if isfield(stats, 'cosine') && isfield(stats.cosine, 'R2')
            R2_cosine = stats.cosine.R2;
        end

        p_perm = NaN;
        if isfield(stats, 'permutation') && isfield(stats.permutation, 'p')
            p_perm = stats.permutation.p;
        end

        p_anova_hue = NaN;
        p_anova_sat = NaN;
        p_anova_int = NaN;
        if isfield(stats, 'anova') && ~isempty(stats.anova)
            if isfield(stats.anova,'p_hue'), p_anova_hue = stats.anova.p_hue; end
            if isfield(stats.anova,'p_sat'), p_anova_sat = stats.anova.p_sat; end
            if isfield(stats.anova,'p_int'), p_anova_int = stats.anova.p_int; end
        end

        peakClass = "";
        if isfield(stats, 'peak') && isfield(stats.peak, 'class')
            peakClass = string(stats.peak.class);
        end

        row = [
            string(U.dateStr), ...
            string(U.unitType), ...
            U.unitID, ...
            U.nTrials, ...
            meanRateOverall, ...
            R_rayleigh, ...
            p_rayleigh, ...
            R2_cosine, ...
            p_perm, ...
            p_anova_hue, ...
            p_anova_sat, ...
            p_anova_int, ...
            peakClass ...
        ];

        fmt = '%s,%s,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%s\n';
        fprintf(fid, fmt, ...
            char(row(1)), char(row(2)), ...
            row(3), row(4), row(5), row(6), row(7), row(8), ...
            row(9), row(10), row(11), row(12), char(row(13)));

        fclose(fid);

    catch ME
        warning('Quick unit CSV summary failed for unit %s_%d: %s', ...
            unitType, unitID, ME.message);
    end
end
