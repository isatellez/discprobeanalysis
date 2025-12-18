function run_vss_spike_accounting(DATES)
% VSS harmonic figures for sat=1 units

clc;

rootDir = fileparts(mfilename('fullpath'));
addpath(genpath(rootDir));

config = load_config();
COL    = load_color_metadata(config);

% handle DATES / "all" style
if nargin < 1
    DATES = [];
end

if ischar(DATES) || isstring(DATES)
    dstr = strtrim(char(DATES));
    if strcmpi(dstr, 'all')
        DATES = [];
    else
        DATES = {dstr};
    end
end

if isfield(config, 'paths') && isfield(config.paths, 'output') && ~isempty(config.paths.output)
    outRoot = config.paths.output;
else
    outRoot = fullfile(rootDir, 'output');
end

vssRoot = fullfile(outRoot, 'vss');

% if no DATES provided, infer all 6-digit date folders under vssRoot
if isempty(DATES)
    if ~isfolder(vssRoot)
        warning('No VSS root folder found at %s', vssRoot);
        return;
    end

    D = dir(vssRoot);
    isDir = [D.isdir];
    names = {D.name};
    names = names(isDir & ~ismember(names, {'.','..'}));

    isDate = cellfun(@(s) ~isempty(regexp(s, '^\d{6}$', 'once')), names);
    dateNames = names(isDate);

    if isempty(dateNames)
        warning('No 6-digit VSS date folders found in %s', vssRoot);
        return;
    end

    DATES = sort(dateNames(:).');
end

DATES = normalize_dates(DATES);

for d = 1:numel(DATES)
    dateStr = char(DATES{d});
    fprintf('\nVSS harmonic figures: %s\n', dateStr);

    sessionRoot = fullfile(vssRoot, dateStr);
    unitsRoot   = fullfile(sessionRoot, 'units');
    figRoot     = fullfile(sessionRoot, 'spikeAccounting');

    if ~isfolder(unitsRoot)
        warning('No VSS units folder found for %s at %s', dateStr, unitsRoot);
        continue;
    end
    if ~isfolder(figRoot)
        mkdir(figRoot);
    end

    Dunits = dir(unitsRoot);
    isDir = [Dunits.isdir];
    names = {Dunits.name};
    keep = isDir & ~ismember(names, {'.','..'});
    Dunits = Dunits(keep);

    for i = 1:numel(Dunits)
        unitLabel = Dunits(i).name;           % e.g. 'SU13', 'MU269'
        unitDir   = fullfile(unitsRoot, unitLabel);

        mats = dir(fullfile(unitDir, 'VSSUnit_*.mat'));
        if isempty(mats)
            continue;
        end

        matPath = fullfile(unitDir, mats(1).name);
        S = load(matPath);

        if ~isfield(S, 'U')
            warning('No U struct in %s, skipping', matPath);
            continue;
        end

        U = S.U;
        if ~isfield(U,'nTrials') || U.nTrials == 0
            continue;
        end

        if isfield(S,'hm')
            hm = S.hm;
        else
            hm = compute_hue_means(U);
        end

        if isfield(S,'P')
            P = S.P;
        else
            P = peak_model(U, config);
        end

        fprintf('  %s: VSS harmonics for %s%d (%d trials, class=%s)\n', ...
            dateStr, U.unitType, U.unitID, U.nTrials, string(P.class));

        % robust plotting: never let a bad fig kill the "all" run
        try
            plot_vss_unit_harmonics(U, hm, P, COL, config, figRoot);
        catch ME
            warning('run_vss_spike_accounting:PlotFailed', ...
                'Failed to plot %s%d on %s: %s', ...
                U.unitType, U.unitID, dateStr, ME.message);
            try
                if ishghandle(gcf)
                    close(gcf);
                end
            catch
            end
        end
    end
end

end
