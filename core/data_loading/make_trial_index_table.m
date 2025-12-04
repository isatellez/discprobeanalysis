function T = make_trial_index_table(S, COL, config)
% Build trial-level table for one session.
% Uses DiscprobeColor and master ProbeCols/ProbeColIDs.

ExptTrialsDisc = S.ExptTrialsDisc;
nTrials        = size(ExptTrialsDisc, 1);
dateStr        = S.dateStr;

if nTrials == 0
    error('make_trial_index_table: no trials found for %s.', dateStr);
end

trial = (1:nTrials).';

% grab trial RGB from DiscprobeColor (first column of ExptTrialsDisc)
try
    trlCols = cell2mat(cellfun(@(s) s.DiscprobeColor(:).', ...
                               ExptTrialsDisc(:,1), ...
                               'UniformOutput', false));
catch
    error('make_trial_index_table: could not read DiscprobeColor from ExptTrialsDisc(:,1) for %s.', dateStr);
end

trlCols   = double(trlCols);

% master probe RGB and IDs from COL
if ~isfield(COL, 'probeCols') || ~isfield(COL, 'probeIDs')
    error('make_trial_index_table: COL.probeCols / COL.probeIDs missing.');
end

masterRGB = double(COL.probeCols(:,1:3));
masterIDs = double(COL.probeIDs);

% direct row-wise matches
[isHit, loc] = ismember(trlCols, masterRGB, 'rows');

trlColIDs = nan(size(trlCols,1), size(masterIDs,2));
trlColIDs(isHit, :) = masterIDs(loc(isHit), :);

% handle non-exact matches with nearest neighbor snapping
need  = find(~isHit);
nnIdx = [];
dists = [];

if ~isempty(need)
    if exist('knnsearch','file')
        [nnIdx, dists] = knnsearch(masterRGB, trlCols(need,:));
    else
        nnIdx = zeros(numel(need),1);
        dists = zeros(numel(need),1);
        for ii = 1:numel(need)
            x   = trlCols(need(ii),:);
            dif = masterRGB - x;
            dst = sqrt(sum(dif.^2,2));
            [dists(ii), nnIdx(ii)] = min(dst);
        end
    end

    % tolerance for snapping in RGB space (tweak if needed)
    tol = 5;
    use = dists <= tol;

    trlColIDs(need(use), :) = masterIDs(nnIdx(use), :);
else
    use   = [];
    dists = [];
end

% bookkeeping for match type and mapping into master list
match_type = strings(nTrials,1);
match_type(:) = "unmatched";

master_row = nan(nTrials,1);
rgb_dist   = nan(nTrials,1);

master_row(isHit) = loc(isHit);
match_type(isHit) = "exact";

if ~isempty(need)
    got = need(use);
    master_row(got) = nnIdx(use);
    match_type(got) = "snap";
    rgb_dist(got)   = dists(use);
end

% unpack IDs and RGB per trial
IDs = trlColIDs;
RGB = trlCols;

% DKL coordinates via spatial_transforms util
dklL = nan(nTrials,1);
dklM = nan(nTrials,1);
dklS = nan(nTrials,1);

if isfield(COL, 'dkl') && ~isempty(COL.dkl)
    [dklL, dklM, dklS] = spatial_transforms.index_to_dkl(master_row, COL);
end

% build trial table
T = table();

T.dateStr = repmat(string(dateStr), nTrials, 1);
T.unitID  = nan(nTrials,1);      % filled in prepare_unit

T.trial   = trial;

T.hueID   = IDs(:,1);            % hue index (1–nHue)
T.satID   = IDs(:,2);            % saturation ID/value
T.lumID   = IDs(:,3);            % luminance/elevation ID

T.dklL    = dklL;
T.dklM    = dklM;
T.dklS    = dklS;

T.matchType = match_type;
T.rgbDist   = rgb_dist;
T.masterRow = master_row;

T.R = RGB(:,1);
T.G = RGB(:,2);
T.B = RGB(:,3);

if isfield(S, 'trTypes')
    T.trType = S.trTypes(:);
else
    T.trType = trial;
end

% optional: save per-session trial index CSV
if isfield(config, 'paths') && isfield(config.paths, 'output') && ~isempty(config.paths.output)
    dateDir = fullfile(config.paths.output, dateStr);
    if ~isfolder(dateDir)
        mkdir(dateDir);
    end
    outCSV = fullfile(dateDir, sprintf('%s_TrialIndex.csv', dateStr));
    try
        writetable(T, outCSV);
    catch ME
        warning('make_trial_index_table: could not write trial index CSV to %s (%s).', outCSV, ME.message);
    end
end

end
