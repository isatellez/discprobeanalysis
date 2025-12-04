function save_unit_outputs(dateDataRoot, U, stats)
% Save per-unit data + stats into a .mat file in a unit-specific folder.

if nargin < 3
    error('save_unit_outputs expects (dateDataRoot, U, stats).');
end

if ~isfield(U, 'unitID') || ~isfield(U, 'unitType') || ~isfield(U, 'dateStr')
    error('U is missing unitID, unitType or dateStr.');
end

if ~exist(dateDataRoot,'dir')
    mkdir(dateDataRoot);
end

% mirror the fig folder naming: SU_123, MU_5, etc.
unitLabel = sprintf('%s_%d', U.unitType, U.unitID);
unitDir   = fullfile(dateDataRoot, unitLabel);
if ~exist(unitDir,'dir')
    mkdir(unitDir);
end

% file tag is consistent with figs
fileTag = sprintf('%s_%s_%d', string(U.dateStr), U.unitType, U.unitID);

matFile = fullfile(unitDir, sprintf('DiscProbeUnit_%s.mat', fileTag));

% light metadata for future sanity checks
meta = struct();
meta.dateStr   = string(U.dateStr);
meta.unitID    = U.unitID;
meta.unitType  = U.unitType;
meta.nTrials   = U.nTrials;
meta.winEarly  = U.winEarly;
meta.winLate   = U.winLate;
meta.winFull   = U.winFull;

try
    save(matFile, 'U', 'stats', 'meta', '-v7');
catch ME
    warning('Could not save unit outputs to %s: %s', matFile, ME.message);
end

% optional: quick CSV row per unit (append to per-date summary within data dir)
try
    csvFile = fullfile(dateDataRoot, sprintf('%s_unit_quickSummary.csv', string(U.dateStr)));
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

    % string -> char, numeric stay numeric; write as CSV line
    fmt = '%s,%s,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%s\n';
    fprintf(fid, fmt, ...
        char(row(1)), char(row(2)), ...
        row(3), row(4), row(5), row(6), row(7), row(8), ...
        row(9), row(10), row(11), row(12), char(row(13)));

    fclose(fid);

catch ME
    warning('Quick unit CSV summary failed for unit %s_%d: %s', ...
        U.unitType, U.unitID, ME.message);
end

end
