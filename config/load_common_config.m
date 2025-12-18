function config = load_common_config()
% load_common_config
% shared settings for all wrappers (discprobe, sphere, wachtler)

% --- global identity ---
% config.monkey = 'Jacomo'; 

% --- paths ---
config.paths.base      = '/Users/tellezi2/Documents/DiscProbe';
config.paths.code      = '/Users/tellezi2/Documents/DiscProbeAnalysis';
config.paths.output    = config.paths.base;
config.paths.global    = '/Users/tellezi2/Documents/DiscProbe';

% --- experimental parameters ---
% spike counting windows
config.time.winEarly   = [0.00 0.20];
config.time.winLate    = [0.25 0.40];
config.time.fullWin    = [0.00 0.40];
config.time.winForHz   = 0.40;

% color space definitions
config.space.mode      = 'dkl';
config.space.zeroHue   = 16; 

% --- data processing defaults ---
config.trials.forceRebuildIndex = false;

% outlier detection 
config.qc.removeOutliers = true; 

% --- plotting defaults ---
config.plot.dpi        = 300;
config.plot.makePlots  = false; 
end