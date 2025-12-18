function build_all_keep_from_wachtler(wachtlerRoot, outKeepPath)
% Build a global ALL_keep.csv from per-session Wachtler fits.
%
% Scans:
%   <DiscProbe>/wachtler/<date>/tables/<date>_wachtler_fits.csv
% and creates:
%   <DiscProbe>/ALL_keep.csv
%
% Unit tags are always SU/MU + phy number (unitID), never unitNum.

config = load_config();

% where the wachtler/<date>/tables live
if nargin < 1 || isempty(wachtlerRoot)
    wachtlerRoot = fullfile(config.paths.base, 'wachtler');
end

% where to write ALL_keep.csv
if nargin < 2 || isempty(outKeepPath)
    outKeepPath = fullfile(config.paths.base, 'ALL_keep.csv');
end

if ~exist(wachtlerRoot, 'dir')
    error('Wachtler root folder not found: %s', wachtlerRoot);
end

% all per-session fits like /wachtler/<date>/tables/<date>_wachtler_fits.csv
D = dir(fullfile(wachtlerRoot, '*', 'tables', '*_wachtler_fits.csv'));
if isempty(D)
    error('No wachtler_fits.csv files found under %s', wachtlerRoot);
end

fprintf('Found %d Wachtler fit files.\n', numel(D));

T = table();

for i = 1:numel(D)
    csvPath = fullfile(D(i).folder, D(i).name);
    [~, dateStr] = fileparts(fileparts(D(i).folder)); % folder above "tables"
    fprintf('  Loading %s...\n', dateStr);

    Ti = readtable(csvPath);
    Ti.session = repmat(string(dateStr), height(Ti), 1);

    % attach unitType (SU/MU) if we can, using per-session unit_summary
    unitTypeCol = repmat("SU", height(Ti), 1); % default if we can't find better

    % first: if fits already have a unitType-ish column, trust it
    vTi = string(Ti.Properties.VariableNames);
    if any(vTi == "unitType")
        unitTypeCol = string(Ti.unitType);
    elseif any(vTi == "unit_type")
        unitTypeCol = string(Ti.unit_type);
    else
        % otherwise, try to read <output>/<date>/tables/<date>_unit_summary.csv
        if isfield(config, 'paths') && isfield(config.paths, 'output') ...
                && ~isempty(config.paths.output)

            sumFile = fullfile(config.paths.output, dateStr, 'tables', ...
                               sprintf('%s_unit_summary.csv', dateStr));
            if isfile(sumFile)
                try
                    Ts = readtable(sumFile);
                    vS = string(Ts.Properties.VariableNames);

                    if any(vS == "unitID") && any(vS == "unitType") ...
                            && any(vTi == "unitID")

                        idSum   = Ts.unitID;
                        typeSum = string(Ts.unitType);
                        idFit   = Ti.unitID;

                        unitTypeCol = repmat("SU", height(Ti), 1);
                        for k = 1:height(Ti)
                            match = (idSum == idFit(k));
                            if any(match)
                                unitTypeCol(k) = typeSum(find(match, 1, 'first'));
                            end
                        end
                    end
                catch ME
                    warning('Could not read unit_summary for %s: %s', dateStr, ME.message);
                end
            end
        end
    end

    Ti.unitType = unitTypeCol;
    T = [T; Ti]; %#ok<AGROW>
end

vn = string(T.Properties.VariableNames);

% session
sess = string(T.session);

% unitType as string
if any(vn == "unitType")
    unitType = string(T.unitType);
elseif any(vn == "unit_type")
    unitType = string(T.unit_type);
else
    unitType = repmat("SU", height(T), 1);
end

% phy number: always from unitID (or phy_id if that ever shows up)
if any(vn == "unitID")
    ids = T.unitID;
elseif any(vn == "phy_id")
    ids = T.phy_id;
else
    error('No unitID or phy_id column found in Wachtler fits tables.');
end

% build SU/MU + phy number tag, e.g. SU13, MU45
unitTag = strings(height(T), 1);
for j = 1:height(T)
    unitTag(j) = sprintf('%s%d', unitType(j), ids(j));
end

% define keep based on available metrics
keep = zeros(height(T), 1);

if any(vn == "isTunedR2")
    keep = T.isTunedR2 ~= 0;
elseif any(vn == "pChi2")
    keep = T.pChi2 < 0.05;
elseif any(vn == "R2")
    keep = T.R2 >= 0.2;
else
    error('No isTunedR2, pChi2, or R2 column found in the fits files.');
end

Tbl = table(sess, unitTag, keep, ...
    'VariableNames', {'session','unit','keep'});

Tbl = unique(Tbl, 'rows');
Tbl = sortrows(Tbl, {'session','unit'});

writetable(Tbl, outKeepPath);
fprintf('\n✅ Wrote ALL_keep.csv with %d units to %s\n', height(Tbl), outKeepPath);
end
