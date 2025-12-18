function config = load_config()

config.sphere.defaultSession= 'Jacomo_250513';


% what figures do we want to do?
% skip tests/figures by setting false
% !!! if you want to drop outliers you have to set these two to true
config.do.boxPlots         = true;
config.do.dropOutliers     = true;

%discprobe "legacy" (the original set of analyses i did)
config.do.tuningFigure     = true; %felix's fig
config.do.rasterCanonical  = true; %raster in rainbow
config.do.rasterByMean     = true; %raster by mean
config.do.permutationTest  = true; %cv permutation test
config.do.spikeAccounting  = true; %cosine fit
config.do.rayleigh         = true; %raylegih
config.do.fourier          = true; % fft
config.do.fourierFig       = true; % to make the 2x2 diagnostic fig
config.do.dklSuite         = true; %polar plots and dome
config.do.swatches         = true; %pantenes

%vss
config.do.peakModel        = true; %fit first and second harmonic

%stas
config.do.stas             = false;

config.plot.makePlots = false;


% Fourier params
config.fourier.maxHarmonics = 8;
config.fourier.useHann      = false;
config.fourier.detrend      = true;

% spike window
config.time.winEarly  = [0.00 0.20];
config.time.winLate   = [0.25 0.40];
config.time.fullWin   = [0.00 0.40];
config.time.winForHz  = 0.40;

% saturation IDs 
config.saturations = [0.33 0.66 1];

% output resolution
config.fig.pngDpi = 300;

% paths 
config.paths.base      = '/Users/tellezi2/Documents/DiscProbe';
config.paths.code      = '/Users/tellezi2/Documents/DiscProbeAnalysis';
config.paths.output    = config.paths.base;


config.trials.forceRebuildIndex = true;


%name ( for filename matching, it needs to be the same) 
% config.monkey = 'Jacomo'; 

%alpha harmonic
config.tuning.alphaHarmonic= 0.05;

% DKL mode or index mode 
config.space.mode = 'dkl';     % or 'index'
config.space.zeroHue = 16;     % which hue sits at 0°

% permutation test settings 
config.pt.nPerm = 1000;
config.pt.alpha = 0.05;
config.pt.seed  = 2025;

% legacy / DiscProbe-v1 style settings
config.legacy.useLegacySats = true;
config.legacy.sats          = [0.33 0.66 1];   % the three main ones
config.legacy.satTol        = 1e-3;            % tolerance for matching

%vss
config.do.saveVssUnitMats  = true;   % turn off per-unit MATs if you want
config.do.buildVssKeep     = true;   % let build_all_keep_vss own ALL_keep_vss.csv
end
