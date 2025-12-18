function config = load_sphere_config()
% config specifically for sphere_slicer.m

% inherit common settings
config = load_common_config();

% --- slicer specifics ---
config.sphere.defaultSession = 'Jacomo_250513';
config.trials.minTrialsSlice = 10;

% --- plotting preferences ---
config.plot.makePlots = false; 
config.plot.suppressSubfolders = true; 
config.plot.saveFigs = false;

% --- analysis flags ---
config.sphere.do.hueMeans        = true;
config.sphere.do.peakModel       = true;
config.sphere.do.cosineFit       = true;
config.sphere.do.rayleighStats   = true;
config.sphere.do.rayleighPlot    = true;
config.sphere.do.spikeAccounting = true;
config.sphere.do.dklSuite = true;

% --- qc ---
config.qc.removeOutliers = true;
end