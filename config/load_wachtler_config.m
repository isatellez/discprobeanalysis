function config = load_wachtler_config()
% config specifically for run_wachtler_analysis.m

% inherit common settings
config = load_common_config();

% --- wachtler specifics ---
config.wachtler.fitGaussian = true;
config.wachtler.fitDog      = false; 

% specific saturations for fitting
config.wachtler.satsToFit   = [0.33 0.66 1.0]; 
end