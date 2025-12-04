function config = load_config()

% what figures do we want to do?
% skip tests/figures by setting false
% !!! if you want to drop outliers you have to set these two to true
config.do.boxPlots         = true;
config.do.dropOutliers     = true;


config.do.tuningFigure     = true;
config.do.rasterCanonical  = true;
config.do.rasterByMean     = true;
config.do.permutationTest  = true;
config.do.spikeAccounting  = true;
config.do.rayleigh         = true;
config.do.fourier          = true;
config.do.fourierQuickFig  = true;
config.do.dklSuite         = true;
config.do.swatches         = true;
config.do.peakModel        = true;
config.do.stas             = true;

% Fourier params
config.fourier.maxHarmonics = 8;
config.fourier.useHann      = true;
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
config.paths.code      = '/Users/tellezi2/Documents/DiscProbesCode';
config.paths.output    = config.paths.base;

% monkey name ( for filename matching, it needs to be the same) 
config.monkey = 'Jacomo'; %jocamo has been renamed

% DKL mode or index mode 
config.space.mode = 'dkl';     % or 'index'
config.space.zeroHue = 16;     % which hue sits at 0°

% permutation test settings 
config.pt.nPerm = 1000;
config.pt.alpha = 0.05;
config.pt.seed  = 2025;

end
