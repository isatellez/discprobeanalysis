function run_fourier_population_analysis(MODE, dateToAnalyze, spaceTag)
% Fourier population keeper + analysis (single / pooled / VSS)

DO_PANEL_A        = true;
DO_PANEL_B        = true;
DO_UNIT_OVERLAYS  = false;
DO_POWER_SPECTRUM = true;

if nargin < 1 || isempty(MODE)
    MODE = 'all';                          % 'single', 'all', or 'vss'
end
if nargin < 2
    dateToAnalyze = '';
end
if nargin < 3 || isempty(spaceTag)
    spaceTag = 'DKL';
end

clc;
set(0,'DefaultFigureVisible','off');

% make sure repo is on path
rootDir = fileparts(mfilename('fullpath'));
addpath(genpath(rootDir));

config   = load_config();
basePath = config.paths.base;             % e.g. /Users/.../DiscProbe
rgbPath  = config.paths.code;             % where rgb.csv lives

fourSubdir      = 'overall';
fourBase        = fullfile(basePath, 'Fourier');
fourFilePattern = sprintf('%%s_FourierHue_summary_%s.csv', spaceTag);

if ~exist(fourBase,'dir')
    mkdir(fourBase);
end

if strcmpi(MODE,'single') && (isempty(dateToAnalyze))
    error('For MODE=''single'', you need to pass dateToAnalyze, e.g. ''250513''.');
end

rgbCSV = fullfile(rgbPath, 'rgb.csv');
pal16  = load_palette16(rgbCSV);
if isempty(pal16)
    warning('rgb.csv missing or unreadable -> skipping hue ring.');
else
    fprintf('palette OK: min=%.3f max=%.3f\n', min(pal16(:)), max(pal16(:)));
end

% I/O root
fourBase = fullfile(basePath,'Fourier');
if ~exist(fourBase,'dir'), mkdir(fourBase); end

%% ========================== choose units (MODE) ==========================

switch lower(MODE)
    case 'single'
        fourFile = sprintf(fourFilePattern, dateToAnalyze);
        fourCSV  = fullfile(basePath, dateToAnalyze, fourSubdir, fourFile);
        if ~isfile(fourCSV)
            error('Missing Fourier file: %s', fourCSV);
        end

        outDir   = fullfile(fourBase, dateToAnalyze);
        if ~exist(outDir,'dir'), mkdir(outDir); end
        keepPath = fullfile(outDir, sprintf('%s_keep.csv', dateToAnalyze));

        T = readtable(fourCSV);
        T.session = repmat(string(dateToAnalyze), height(T), 1);

        [ordType, ordNum, tag] = unit_sort_keys(T);
        T.unit_tag  = tag; 
        T.typeOrder = ordType; 
        T.idNum     = ordNum;
        T           = sortrows(T, {'session','typeOrder','idNum','unit_tag'});

        if ~isfile(keepPath)
            Tkeep = make_keep_table_from_unit(T);
            writetable(Tkeep, keepPath);
            fprintf('Created keeper file: %s\nEdit keep=1 then rerun.\n', keepPath);
            return
        else
            Tkeep = readtable(keepPath);
            T     = innerjoin(T, Tkeep, 'Keys', {'session','unit'});
            Tuse  = T(T.keep==1, :);
            if isempty(Tuse)
                warning('No units marked keep=1. Nothing to analyze.');
                return
            end
            [ordType, ordNum, tag] = unit_sort_keys(Tuse);
            Tuse.unit_tag  = tag; 
            Tuse.typeOrder = ordType; 
            Tuse.idNum     = ordNum;
            Tuse = sortrows(Tuse, {'session','typeOrder','idNum','unit_tag'});
        end

    case 'all'
        % pooled analysis, keeper file in Fourier/ALL/ALL_keep.csv
        searchPattern = sprintf('*_FourierHue_summary_%s.csv', spaceTag);
        D = dir(fullfile(basePath, '*', fourSubdir, searchPattern));
        if isempty(D)
            error('No FourierHue_summary_%s.csv files found under %s/*/%s/', ...
                  spaceTag, basePath, fourSubdir);
        end

        sess  = strings(numel(D),1);
        paths = strings(numel(D),1);
        for i = 1:numel(D)
            paths(i) = string(fullfile(D(i).folder, D(i).name));
            parentOfSubdir = fileparts(D(i).folder);      % /basePath/<session>
            [~, sess(i)]   = fileparts(parentOfSubdir);   % <session> from that folder name
        end

        outDir   = fullfile(fourBase, 'ALL');
        if ~exist(outDir,'dir'), mkdir(outDir); end
        keepPath = fullfile(outDir, 'ALL_keep.csv');

        T = harmonize_fourier_table(paths, sess);

        [ordType, ordNum, tag] = unit_sort_keys(T);
        T.unit_tag  = tag; 
        T.typeOrder = ordType; 
        T.idNum     = ordNum;
        T           = sortrows(T, {'session','typeOrder','idNum','unit_tag'});

        if ~isfile(keepPath)
            Tkeep = make_keep_table_from_unit(T);
            writetable(Tkeep, keepPath);
            fprintf('Created pooled keeper: %s\nEdit keep=1 then rerun.\n', keepPath);
            return
        else
            Tkeep = readtable(keepPath);
            T     = innerjoin(T, Tkeep, 'Keys', {'session','unit'});
            Tuse  = T(T.keep==1, :);
            if isempty(Tuse)
                warning('No units marked keep=1. Nothing to analyze.');
                return
            end
            [ordType, ordNum, tag] = unit_sort_keys(Tuse);
            Tuse.unit_tag  = tag; 
            Tuse.typeOrder = ordType; 
            Tuse.idNum     = ordNum;
            Tuse = sortrows(Tuse, {'session','typeOrder','idNum','unit_tag'});
        end

        case 'vss'
        % pooled analysis using VSS Fourier CSVs + ALL_keep_vss.csv

        % where VSS outputs live
        if isfield(config,'paths') && isfield(config.paths,'output') && ~isempty(config.paths.output)
            vssRoot = fullfile(config.paths.output, 'vss');
        else
            vssRoot = fullfile(basePath, 'vss');
        end

        % find all per-date VSS Fourier summaries
        searchPattern = '*_vss_fourier_sat1.csv';
        D = dir(fullfile(vssRoot, '*', 'tables', searchPattern));
        if isempty(D)
            error('No VSS Fourier files (%s) found under %s/*/tables/', ...
                  searchPattern, vssRoot);
        end

        % where population outputs go
        outDir = fullfile(fourBase, 'VSS');
        if ~exist(outDir,'dir'), mkdir(outDir); end

        % VSS keep list (built by run_vss_peakmodels)
        keepPath = fullfile(basePath, 'ALL_keep_vss.csv');
        if ~isfile(keepPath)
            error('VSS keep list not found at %s. Run run_vss_peakmodels(''all'') first.', keepPath);
        end

        % build pooled Fourier table T from VSS CSVs
        T = table();
        for i = 1:numel(D)
            csvPath = fullfile(D(i).folder, D(i).name);

            % session folder is parent of "tables"
            sessFolder = fileparts(D(i).folder);   % .../vss/<date>
            [~, sess]  = fileparts(sessFolder);    % <date> (e.g. 241122)

            Ti = readtable(csvPath);
            Ti.session = repmat(string(sess), height(Ti), 1);

            % construct a "unit" column that matches ALL_keep_vss.unit
            % prefer unitLabel if present; else build from type + ID
            vn = string(Ti.Properties.VariableNames);
            if any(vn == "unitLabel")
                unitCol = string(Ti.unitLabel);
            else
                uType = strings(height(Ti),1);
                if any(vn == "unitType")
                    uType = string(Ti.unitType);
                end
                uID = strings(height(Ti),1);
                if any(vn == "unitID")
                    uID = string(Ti.unitID);
                end
                unitCol = uType + uID;   % "SU13", "MU4", etc.
            end
            Ti.unit = unitCol;

            % append
            T = [T; Ti]; %#ok<AGROW>
        end

        % sort units in a nice order
        [ordType, ordNum, tag] = unit_sort_keys(T);
        T.unit_tag  = tag;
        T.typeOrder = ordType;
        T.idNum     = ordNum;
        T           = sortrows(T, {'session','typeOrder','idNum','unit_tag'});

        % apply ALL_keep_vss.csv
        Tkeep = readtable(keepPath);
        if any(strcmpi(Tkeep.Properties.VariableNames, 'keep'))
            Tkeep = Tkeep(Tkeep.keep == 1, :);
        end

        % intersection: only units present in BOTH the VSS Fourier table and keep list
        Tuse = innerjoin(T, Tkeep, 'Keys', {'session','unit'});
        if isempty(Tuse)
            warning('No overlap between VSS Fourier table and ALL_keep_vss keep==1. Nothing to analyze.');
            return
        end

        % re-sort after join
        [ordType, ordNum, tag] = unit_sort_keys(Tuse);
        Tuse.unit_tag  = tag;
        Tuse.typeOrder = ordType;
        Tuse.idNum     = ordNum;
        Tuse = sortrows(Tuse, {'session','typeOrder','idNum','unit_tag'});

    otherwise
        error('MODE must be ''single'', ''all'', or ''vss''.');
end

%% ========================== core analysis (shared) =========================

v = string(Tuse.Properties.VariableNames);
pick = @(alts) alts(find(ismember(alts, v), 1, 'first'));

col_amp1   = pick(["amp1","over_amp1","k1_amp"]);
col_amp2   = pick(["amp2","over_amp2","k2_amp"]);
col_phase1 = pick(["phase1","over_phase1","k1_phase","phase1_deg"]);
col_phase2 = pick(["phase2","over_phase2","k2_phase","phase2_deg"]);
col_DC     = pick(["DC","over_DC","dc"]);
if any(cellfun(@isempty, {col_amp1,col_amp2,col_phase1,col_phase2,col_DC}))
    error('Table must contain k1/k2 amplitude & phase and DC (tried name variants).');
end

amp1   = Tuse.(col_amp1);
amp2   = Tuse.(col_amp2);
ph1    = mod(Tuse.(col_phase1), 360);
ph2raw = mod(Tuse.(col_phase2), 360);
ph2    = mod(ph2raw, 180);
DCv    = max(Tuse.(col_DC), eps);

% --- override with true preferred angles if available (VSS pipeline) ---
if any(v == "pref_deg_h1")
    ph1 = mod(Tuse.pref_deg_h1, 360);
end

if any(v == "axis_deg_h2")
    % axis is already 0–180; for axial plots we want orientation, not raw phase
    ph2raw = mod(Tuse.axis_deg_h2, 360);  % for axial histogram
    ph2    = mod(Tuse.axis_deg_h2, 180);  % for width vs orientation
end


% widths from each unit's mean-per-hue curve
nUnits = height(Tuse);
width_long  = nan(nUnits,1);
width_short = nan(nUnits,1);
for i = 1:nUnits
    sess_i = string(Tuse.session(i));
    tag_i  = derive_unit_tag(Tuse, i);
    [ok, ybar] = get_unit_ybar(basePath, sess_i, tag_i);
    if ~ok || isempty(ybar), continue; end
    [w_long, w_short] = fwhm_circular_from_ybar(ybar);
    width_long(i)  = w_long;
    width_short(i) = w_short;
end

% dominance split
idx_uni = amp1 > amp2;
idx_bi  = amp2 > amp1;

% saver tag
if strcmpi(MODE,'single')
    tagRun = string(dateToAnalyze);
elseif strcmpi(MODE,'vss')
    tagRun = "VSS";
else
    tagRun = "ALL";
end
sav = @(fn) fullfile(outDir, sprintf('%s_%s', tagRun, fn));

% ==== unified angular bin definition (shared by hist + hue ring) ====
nBins = 16;
d     = 360 / nBins;
binEdges_deg = (-d/2) : d : (360 - d/2);

opts = struct('offset_deg', 0, 'flip', false);

% phase vs width
fW = figure('Color','w','Position',[80 80 1000 430]);
subplot(1,2,1);
ii = find(idx_uni & ~isnan(width_long));
scatter(ph1(ii), width_long(ii), 46, 'k','filled'); grid on
xlabel('Phase k=1 (deg)'); ylabel('FWHM width (deg)'); title('k=1 dominant');
subplot(1,2,2);
jj = find(idx_bi & ~isnan(width_short));
scatter(ph2(jj), width_short(jj), 46, 'k','filled'); grid on
xlabel('Phase k=2 (deg, 0–180)'); ylabel('FWHM width (deg)'); title('k=2 dominant');
saveas(fW, sav('PhaseVsWidth_k1_k2.png')); close(fW);

% Panel A
if DO_PANEL_A
    fA = figure('Color','w','Position',[400 300 1800 900]);

    ax1 = polaraxes('Parent', fA);
    ax1.Position = [0.14 0.12 0.32 0.76];

    polarhist_counts_with_mean(ph1(idx_uni), binEdges_deg, false, ...
        'First harmonic dominant', ax1, pal16, opts);

    ax2 = polaraxes('Parent', fA);
    ax2.Position = [0.54 0.12 0.32 0.76];

    polarhist_counts_with_mean(ph2raw(idx_bi), binEdges_deg, true, ...
        'Second harmonic dominant', ax2, pal16, opts);

    saveas(fA, sav('PanelA_k1vsK2_dominance.png'));
    close(fA);
end

% Panel B
if DO_PANEL_B
    fB = figure('Color','w','Position',[400 300 1800 900]);

    ax3 = polaraxes('Parent', fB);
    ax3.Position = [0.14 0.12 0.32 0.76];
    polarhist_counts_with_mean(ph1, binEdges_deg, false, ...
        'All cells — 1st harmonic', ax3, pal16, opts);

    ax4 = polaraxes('Parent', fB);
    ax4.Position = [0.54 0.12 0.32 0.76];
    polarhist_counts_with_mean(ph2raw, binEdges_deg, true, ...
        'All cells — 2nd harmonic (axial)', ax4, pal16, opts);

    saveas(fB, sav('PanelB_allCells_k1_and_k2.png'));
    close(fB);
end

% CSV: phases + widths
outCSV = sav('KeptUnits_phase_widths.csv');
S = table( ...
    Tuse.session, Tuse.unit, ...
    amp1, ph1, width_long, ...
    amp2, ph2, width_short, DCv, ...
    'VariableNames', {'session','unit', ...
    'amp1','phase1_deg','width_FWHM_deg', ...
    'amp2','phase2_deg_0to180','width_FWHM_deg_ph2','DC'});
writetable(S, outCSV);
fprintf('Saved analysis to: %s\n', outCSV);

% per-unit overlays
if DO_UNIT_OVERLAYS
    hasUnitFile = ismember('unit_file', string(Tuse.Properties.VariableNames));
    for i = 1:nUnits
        sess_i = string(Tuse.session(i));
        if hasUnitFile && ~isempty(Tuse.unit_file{i})
            tag_i = string(Tuse.unit_file{i});
            if endsWith(tag_i,".mat"), tag_i = erase(tag_i,".mat"); end
        else
            tag_i = derive_unit_tag(Tuse, i);
        end
        try
            plot_fourier_fit_unit_k1k2(basePath, sess_i, tag_i);
        catch ME
            warning('Overlay failed for %s / %s: %s', sess_i, tag_i, ME.message);
        end
    end
end

% population power spectrum (no DC)
if DO_POWER_SPECTRUM
    P_all = [];
    sess_list = strings(0,1);
    unit_list = strings(0,1);
    N0 = NaN;

    for i = 1:nUnits
        sess_i = string(Tuse.session(i));
        tag_i  = derive_unit_tag(Tuse, i);
        [ok, ybar] = get_unit_ybar(basePath, sess_i, tag_i);
        if ~ok || isempty(ybar), continue; end
        y = double(ybar(:));
        if any(~isfinite(y)) || numel(y) < 4, continue; end
        if isnan(N0), N0 = numel(y); end
        if numel(y) ~= N0
            warning('Skipping %s (N mismatch)', tag_i);
            continue;
        end

        y = y - mean(y);
        N = numel(y);
        X = fft(y);
        Kmax = floor(N/2);
        Ak = (2/N) * abs(X(2:Kmax+1));
        Pk = Ak.^2;
        s  = sum(Pk);
        if s <= 0, continue; end

        P_all(end+1,1:Kmax) = (Pk(:)/s).'; %#ok<AGROW>
        sess_list(end+1,1)  = sess_i;      %#ok<AGROW>
        unit_list(end+1,1)  = tag_i;       %#ok<AGROW>
    end

    if isempty(P_all)
        warning('No per-unit power spectra computed. Nothing to plot.');
    else
        meanP  = mean(P_all, 1);
        k      = 1:size(P_all,2);
        nSpec  = size(P_all,1);

        f = figure('Color','w','Position',[120 120 560 420]);
        plot(k, meanP, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'k'); grid on
        xlim([0.5, numel(k)+0.5]); xticks(k);
        xlabel('Harmonic k'); ylabel('Normalized power (no DC)');
        title(sprintf('%s: Population power spectrum (n=%d kept units; %d with valid spectra)', ...
              char(tagRun), nUnits, nSpec));
        saveas(f, sav('PowerSpectrum_noDC.png')); close(f);

        outUnitCSV = sav('PowerSpectrum_perUnit_noDC.csv');
        UnitTab = table(sess_list, unit_list, 'VariableNames', {'session','unit_tag'});
        Ktab = array2table(P_all, 'VariableNames', compose('k%d', k));
        writetable([UnitTab Ktab], outUnitCSV);

        outPopCSV = sav('PowerSpectrum_population_noDC.csv');
        PopTab = table(k.', meanP.', 'VariableNames', {'k','mean_norm_power'});
        writetable(PopTab, outPopCSV);

        fprintf('Saved power spectrum figure + CSVs in %s\n', outDir);
    end
end

%% ======================== helpers (local functions) ========================







% ------------------- rest of helpers unchanged -------------------

function tag_i = derive_unit_tag(T, i)
v = string(T.Properties.VariableNames);
if any(v=="unit_file") && ~isempty(T.unit_file{i})
    tag_i = string(erase(T.unit_file{i}, ".mat"));
elseif any(v=="unit_label") && ~isempty(T.unit_label{i})
    tag_i = string(T.unit_label{i});
else
    u = T.unit(i);
    if iscell(u), u = u{1}; end
    if isstring(u) || ischar(u)
        tag_i = string(u);
    elseif isnumeric(u)
        tag_i = "SU" + string(u);
    else
        tag_i = string(u);
    end
end
tag_i = strtrim(tag_i);
end

function [ok, ybar] = get_unit_ybar(basePath, date6, unitTag)
ok = false; ybar = [];
meanCSV = fullfile(basePath, char(date6), 'Figs', 'MeanNormPerHue_perUnit.csv');
if ~isfile(meanCSV), return; end
T = readtable(meanCSV);
vn = string(T.Properties.VariableNames);
isHueCol = startsWith(vn,"over_h_") | startsWith(vn,"over_hue_") | startsWith(vn,"over_h");
hueCols  = vn(isHueCol);
if isempty(hueCols), return; end
getN = @(s) sscanf(regexprep(s,'[^0-9]',''),'%d');
nums = arrayfun(getN, hueCols); [~, ord] = sort(nums); hueCols = hueCols(ord);

rowMatch = [];
dateStrs = string(T.date);
candDate = dateStrs == string(date6);
m = regexp(char(unitTag),'^(SU|MU)(\d+)$','tokens','once');
if ~isempty(m)
    idNum = str2double(m{2});
    if any(vn=="phy_id")
        cand = candDate & (T.phy_id == idNum);
        if any(cand), rowMatch = find(cand,1,'first'); end
    end
    if isempty(rowMatch) && any(vn=="unit_num")
        cand = candDate & (T.unit_num == idNum);
        if any(cand), rowMatch = find(cand,1,'first'); end
    end
end
if isempty(rowMatch)
    if any(vn=="unit_file")
        uf = erase(string(T.unit_file), ".mat");
        cand = candDate & strcmpi(uf, unitTag);
        if any(cand), rowMatch = find(cand,1,'first'); end
    elseif any(vn=="unit_label")
        cand = candDate & strcmpi(string(T.unit_label), unitTag);
        if any(cand), rowMatch = find(cand,1,'first'); end
    end
end
if isempty(rowMatch)
    cand = find(candDate);
    if ~isempty(cand), rowMatch = cand(1); end
end
if isempty(rowMatch), return; end

ybar = double(T{rowMatch, hueCols}).';
ok = true;
end

function [w_long, w_short] = fwhm_circular_from_ybar(ybar)
N = numel(ybar);
theta = linspace(0, 360, N+1); theta(end) = [];
fineN = 3600;
fineTheta = linspace(0, 360, fineN+1); fineTheta(end) = [];
yf = interp1(theta, ybar(:), fineTheta, 'pchip');
pmax = max(yf);
if ~isfinite(pmax) || pmax<=0
    w_long = NaN; w_short = NaN; return;
end
halfVal = pmax/2;
above = yf >= halfVal;
if ~any(above)
    w_long = 0; w_short = 0; return;
end
ab2 = [above, above];
d  = diff([false, ab2, false]);
starts = find(d==1);
ends   = find(d==-1) - 1;
mask = starts <= fineN;
starts = starts(mask); ends = ends(mask);
lens = ends - starts + 1;
degPerIdx = 360 / fineN;
widths = lens * degPerIdx;
if isempty(widths)
    w_long = 0; w_short = 0;
else
    w_long  = max(widths);
    w_short = min(widths);
end
end

function [typeOrder, idNum, tag] = unit_sort_keys(T)
v = string(T.Properties.VariableNames);
if any(v=="unit_file")
    tag = string(erase(T.unit_file, ".mat"));
elseif any(v=="unit_label")
    tag = string(T.unit_label);
else
    u = T.unit;
    if iscell(u)
        tag = string(u);
    elseif isstring(u) || ischar(u)
        tag = string(u);
    elseif isnumeric(u)
        tag = "SU" + string(u);
    else
        tag = string(u);
    end
end
tag = strtrim(tag);
up  = upper(tag);
isSU = startsWith(up, "SU");
isMU = startsWith(up, "MU");
typeOrder = ones(height(T),1) * 2;
typeOrder(isSU) = 0;
typeOrder(isMU) = 1;
numCell = regexp(tag, '\d+', 'match');
idNum = nan(height(T),1);
for i = 1:height(T)
    if ~isempty(numCell{i})
        idNum(i) = str2double(numCell{i}{1});
    end
end
end

function T = harmonize_fourier_table(paths, sess)
baseNum = ["amp1","over_amp1","k1_amp", ...
           "amp2","over_amp2","k2_amp", ...
           "phase1","over_phase1","k1_phase","phase1_deg", ...
           "phase2","over_phase2","k2_phase","phase2_deg", ...
           "DC","over_DC","dc", ...
           "phy_id","unit_num"];
baseStr = ["session","unit","unit_file","unit_label","date"];
T = table();
for i = 1:numel(paths)
    Ti = readtable(paths(i));
    Ti.session = repmat(sess(i), height(Ti), 1);
    To = table();
    for c = baseStr
        if any(strcmp(Ti.Properties.VariableNames, c))
            To.(c) = string(Ti.(c));
        else
            To.(c) = strings(height(Ti),1);
        end
    end
    for c = baseNum
        if any(strcmp(Ti.Properties.VariableNames, c))
            To.(c) = double(Ti.(c));
        else
            To.(c) = nan(height(Ti),1);
        end
    end
    T = [T; To]; %#ok<AGROW>
end
end

function Tkeep = make_keep_table_from_unit(T)
Tuniq = unique(T(:,{'session','unit'}), 'rows');

u  = string(Tuniq.unit);
up = upper(u);

isSU = startsWith(up, "SU");
isMU = startsWith(up, "MU");

typeOrder = ones(height(Tuniq),1) * 2;
typeOrder(isSU) = 0;
typeOrder(isMU) = 1;

numCell = regexp(u, '\d+', 'match');
idNum   = nan(height(Tuniq),1);
for i = 1:height(Tuniq)
    if ~isempty(numCell{i})
        idNum(i) = str2double(numCell{i}{1});
    end
end

Tuniq.typeOrder = typeOrder;
Tuniq.idNum     = idNum;

Tuniq = sortrows(Tuniq, {'session','typeOrder','idNum','unit'});

Tuniq.keep = zeros(height(Tuniq),1);
Tkeep = removevars(Tuniq, {'typeOrder','idNum'});
end

end

function pal16 = load_palette16(csvPath)
% Load 16 full-saturation RGB colors for the hue ring from rgb.csv.

pal16 = [];

if nargin < 1 || isempty(csvPath)
    warning('load_palette16:NoPath', 'No csvPath provided.');
    return;
end

if ~isfile(csvPath)
    warning('load_palette16:MissingFile', 'rgb.csv not found at %s', csvPath);
    return;
end

A = readmatrix(csvPath);
if isempty(A) || size(A,2) < 3
    warning('load_palette16:BadFile', 'rgb.csv has wrong shape or is empty.');
    return;
end

% assume every 3rd row is full-sat: 3,6,9,...
idx = 3:3:size(A,1);
if isempty(idx)
    warning('load_palette16:NoFullSatRows', ...
        'No full-saturation rows (3:3:end) found in %s', csvPath);
    return;
end

% we want 16 hues if possible
if numel(idx) < 16
    warning('load_palette16:FewHues', ...
        'Only %d full-sat hues found, need 16. Using what is available.', numel(idx));
    idx = idx(1:numel(idx));
else
    idx = idx(1:16);
end

P = double(A(idx,1:3));
P(~isfinite(P)) = 0;

% normalize:
%  - if max > 2, assume 0–255 and divide by 255
%  - if 1 < max <= 2, divide by max (weird 0–1 plus rogue 256 case)
mx = max(P(:));
if mx > 2
    P = P / 255;
elseif mx > 1
    P = P / mx;
end

% clip to [0,1] just in case
P = max(0, min(P, 1));

% rotate so the last hue sits at 0°, i.e., [16 1 2 ... 15]
if size(P,1) >= 16
    P = circshift(P, 1, 1);
end

pal16 = P;
end

function draw_hue_ring(ax, pal16, ring_in, ring_out, edges_deg)
% pal16: 16×3 RGB in [0,1], row 1 = bin 1 color
% ring_in / ring_out: inner/outer radius of the ring
% edges_deg: 1×(N+1) bin edges used for the histogram

if isempty(pal16) || isempty(ax) || ~isvalid(ax)
    return;
end

edges_deg = edges_deg(:).';
N = numel(edges_deg) - 1;

if size(pal16,1) > N
    pal16 = pal16(1:N,:);
end

pal16(~isfinite(pal16)) = 0;
pal16 = max(0, min(pal16,1));

hp = gobjects(N,1);
for k = 1:N
    t1 = deg2rad(edges_deg(k));
    t2 = deg2rad(edges_deg(k+1));
    tt = linspace(t1, t2, 40);

    th_poly = [tt, fliplr(tt)];
    r_poly  = [ring_in * ones(size(tt)), ...
               fliplr(ring_out * ones(size(tt)))];

    hp(k) = patch(ax, th_poly, r_poly, pal16(k,:), ...
        'EdgeColor','none', ...
        'HandleVisibility','off');
end

uistack(hp,'bottom');
end

function polarhist_counts_with_mean(theta_deg, nBins_or_edges, isAxial, ttl, ax, pal16, opts)
% opts.offset_deg (default 0) shifts phase origin
% opts.flip       (default false) flips direction (clockwise vs CCW)

if nargin < 7 || isempty(opts), opts = struct; end
if ~isfield(opts,'offset_deg'), opts.offset_deg = 0; end
if ~isfield(opts,'flip'),       opts.flip       = false; end

sgn = opts.flip * (-1) + (~opts.flip) * 1;

% ---------- handle nBins vs explicit bin edges ----------
if isscalar(nBins_or_edges)
    nBins = nBins_or_edges;
    edges_deg = linspace(0, 360, nBins+1);
else
    edges_deg = nBins_or_edges(:).';
    nBins     = numel(edges_deg) - 1;
end

% ---------- clean angles ----------
theta_deg = theta_deg(:);
theta_deg = theta_deg(isfinite(theta_deg));
nCells    = numel(theta_deg);

theta_deg = mod(theta_deg, 360);
theta_deg = mod(opts.offset_deg + sgn*theta_deg, 360);

hold(ax,'on');

% ---------- directional vs axial ----------
if isAxial
    t = mod(theta_deg, 180);          % 0..180
    theta_plot = [t; t+180];          % duplicate 180° apart
    w = ones(size(theta_plot));
else
    theta_plot = theta_deg;
    w = ones(size(theta_deg));
end

% ---------- binning ----------
binIdx = discretize(theta_plot, edges_deg);
valid  = ~isnan(binIdx);
counts = accumarray(binIdx(valid), w(valid), [nBins 1], @sum, 0).';

% enforce perfect mirroring for axial case
if isAxial
    half = nBins/2;
    for k = 1:half
        m = 0.5 * (counts(k) + counts(k+half));
        counts(k)      = m;
        counts(k+half) = m;
    end
end

% ---------- radial layout ----------
maxC = max(counts);
if maxC <= 0, maxC = 1; end

if maxC <= 5
    step = 1;
elseif maxC <= 10
    step = 2;
else
    step = ceil(maxC/5);
end
outerGrid = step * ceil(maxC/step);

rt = 0:step:outerGrid;
if rt(end) ~= outerGrid
    rt = [rt outerGrid];
end
rticks(ax, rt);
ax.RTickLabel = string(ax.RTick);

ring_gap   = 0.10 * outerGrid;
ring_thick = 0.18 * outerGrid;

ring_in  = outerGrid + ring_gap;
ring_out = ring_in + ring_thick;
outerR   = ring_out + 0.1 * outerGrid;
rlim(ax, [0, outerR]);

% ---------- histogram ----------
polarhistogram(ax, ...
    'BinEdges',  deg2rad(edges_deg), ...
    'BinCounts', counts, ...
    'FaceColor', [0.6 0.6 0.6], ...
    'EdgeColor', 'none');

% ---------- theta grid / spokes ----------
theta_spokes = edges_deg(1:end-1);

ax.ThetaZeroLocation = 'right';          % 0° at +[L−M]
ax.ThetaDir          = 'counterclockwise';
ax.ThetaTick         = theta_spokes;
ax.ThetaTickLabel    = repmat({''}, size(theta_spokes));
ax.ThetaGrid         = 'on';
ax.RGrid             = 'on';
ax.GridAlpha         = 0.25;
ax.Layer             = 'top';

% ---------- numeric labels ----------
majorAngles = [0 90 180 270];
rNum = ring_out + 0.05 * outerGrid;

for k = 1:numel(majorAngles)
    ang = majorAngles(k);
    th  = deg2rad(ang);

    hal = 'center';
    val = 'middle';
    switch ang
        case 90
            val = 'bottom';
        case 270
            val = 'top';
    end

    text(ax, th, rNum, sprintf('%d°', ang), ...
        'HorizontalAlignment', hal, ...
        'VerticalAlignment',   val, ...
        'FontSize',            9);
end

% ---------- DKL axis labels ----------
labelAngles = [0 90 180 270];
labelText   = { '+[L−M]', '+[S−(L+M)]', '−[L−M]', '−[S−(L+M)]' };

rLabel = ring_out + 0.18 * outerGrid;

for k = 1:numel(labelAngles)
    ang = labelAngles(k);
    th  = deg2rad(ang);

    hal = 'center';
    val = 'middle';
    rot = 0;

    switch ang
        case 0
            rot = -90;
        case 180
            rot = 90;
        case 90
            rot = 0;
            val = 'top';
        case 270
            rot = 0;
            val = 'bottom';
    end

    text(ax, th, rLabel, labelText{k}, ...
        'HorizontalAlignment', hal, ...
        'VerticalAlignment',   val, ...
        'FontWeight',          'bold', ...
        'Rotation',            rot);
end

% ---------- title ----------
ht = title(ax, sprintf('%s (n=%d)', ttl, nCells), 'FontWeight', 'bold');
ht.Units = 'normalized';
pos = ht.Position;
pos(2) = 1.06;
ht.Position = pos;

% ---------- hue ring ----------
if nargin >= 6 && ~isempty(pal16)
    draw_hue_ring(ax, pal16, ring_in, ring_out, edges_deg);
end
end
