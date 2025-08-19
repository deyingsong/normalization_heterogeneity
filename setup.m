%% Setup script for EIF network simulation

% Define project root and key directories
rootDir   = pwd;
dataDir   = fullfile(rootDir, 'data');
utilsDir  = fullfile(rootDir, 'src', 'utils');
figDir    = fullfile(rootDir, 'src', 'figure');
simDir    = fullfile(rootDir, 'src', 'simulation');
simScriptDir = fullfile(rootDir, 'sim_script');
analysisDir  = fullfile(utilsDir, '+utils_analysis');

% Add source directories to MATLAB path
addpath(utilsDir);
addpath(figDir);
addpath(simDir);

%% Compile MEX functions

compileFile = fullfile(simScriptDir, 'compile_all.m');
run(compileFile);
mex('-outdir', analysisDir, fullfile(analysisDir, 'spktime2count.c'));

