function config = load_discprobe_config()
% config specifically for run_discprobe_analysis.m

% inherit common settings
config = load_common_config();

% --- analysis flags ---
config.do.tuningFigure     = true;  % felix's fig
config.do.rasterCanonical  = false;  % rainbow raster
config.do.rasterByMean     = false; 
config.do.spikeAccounting  = true;  % cosine fit viz
config.do.rayleigh         = true;  % rayleigh viz
config.do.boxPlots         = true;
config.do.fourierFig       = true;  
config.do.dklSuite         = true;  % polar plots

% --- statistics settings ---
config.do.permutationTest  = true;
config.do.fourier          = true;
config.do.peakModel        = true; 

% --- algorithm parameters ---
% fourier
config.fourier.maxHarmonics = 8;
config.fourier.useHann      = false;
config.fourier.detrend      = true;

% permutation test
config.pt.nPerm = 1000;
config.pt.alpha = 0.05;
config.pt.seed  = 2025;

% tuning fit
config.tuning.alphaHarmonic = 0.05;

% --- legacy support ---
config.legacy.useLegacySats = true;
config.legacy.sats          = [0.33 0.66 1];
end