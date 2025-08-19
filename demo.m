%% DEMO: One-click orientation maps, filters, weights, and a short simulation
% Purpose
%   Let anyone (students to experts) run a quick end-to-end demo:
%   1) Generate V1 & V4/MT orientation maps
%   2) Build demo Gabor filters for V1
%   3) Create demo feedforward & recurrent weights
%   4) Preview V1 rates (toy)
%   5) Run a brief spiking simulation and visualize MT rates + a raster clip
%
% How to use
%   - Open this file in MATLAB and press "Run". No inputs required.
%   - Results (.mat and .pdf/.fig where applicable) are saved under ./data and ./result.
%
% Requirements (typical for this repo)
%   - MATLAB R2021b+ with Statistics, Image Processing, Curve Fitting Toolboxes (if used elsewhere)
%   - A supported C compiler for MEX (see README). Ensure MEX files are built (e.g., via compile_all).
%

%% ------------------------------------------------------------------------
%  0) Housekeeping & paths
%  ------------------------------------------------------------------------
clear; clc;

% Reproducibility for any randomized pieces inside called functions
rng(1234, 'twister');

repoRoot  = pwd;   

% Common folders 
dataDir   = fullfile(repoRoot, 'data');
resultDir = fullfile(repoRoot, 'result');
srcDir    = fullfile(repoRoot, 'src');

% Add likely source paths 
addpath(fullfile(srcDir, 'simulation'));
addpath(fullfile(srcDir, 'utils'));
addpath(fullfile(srcDir, 'figure'));
addpath(fullfile(srcDir, 'simulation', '+spiking_simulation'));  
addpath(fullfile(repoRoot, 'sim_script'));


%% ------------------------------------------------------------------------
%  1) Parameters
%  ------------------------------------------------------------------------
p = neuronConstants();  % central place for physiological & connectivity settings

% Orientation map sizes and filter constants
Nx_V1  = 100;         % V1 orientation map will be 50 x 100
Nx_V4MT = 200;        % V4/MT orientation map will be 100 x 200
Kc     = 10 * 2 * pi; % orientation map spatial frequency

% Basic V1 rate preview settings
theta1 = 180;   % degrees
theta2 = 90;    % degrees
c1     = 0.5;   % contrast 1
c2     = 0.5;   % contrast 2

% Short simulation parameters
simulation_type = 'default';
T               = 2000;  % ms total duration (kept short for a quick demo)

%% ------------------------------------------------------------------------
%  2) Generate orientation maps (V1 and V4/MT) and save
%  ------------------------------------------------------------------------
disp('Generating orientation maps...');

thetaMap_V1 = physiological_properties.generate_ori_map(Nx_V1, Kc);   % 100x100
thetaMap_V1 = thetaMap_V1(1:50, :);
thetaMap_V1_left     = thetaMap_V1(:, 1:50);
thetaMap_V1_right    = thetaMap_V1(:, 51:100);

save(fullfile(dataDir, 'thetaMap_V1_all.mat'),   'thetaMap_V1');
save(fullfile(dataDir, 'thetaMap_V1_left.mat'),  'thetaMap_V1_left');
save(fullfile(dataDir, 'thetaMap_V1_right.mat'), 'thetaMap_V1_right');

thetaMap_V4MT = physiological_properties.generate_ori_map(Nx_V4MT, Kc);   % 200x200
thetaMap_V4MT = thetaMap_V4MT(1:100, :);
thetaMap_V4MT_left    = thetaMap_V4MT(:, 1:100);
thetaMap_V4MT_right   = thetaMap_V4MT(:, 101:200);

save(fullfile(dataDir, 'thetaMap_V4MT_all.mat'),   'thetaMap_V4MT');
save(fullfile(dataDir, 'thetaMap_V4MT_left.mat'),  'thetaMap_V4MT_left');
save(fullfile(dataDir, 'thetaMap_V4MT_right.mat'), 'thetaMap_V4MT_right');

figure;
imagesc(thetaMap_V1);
pbaspect([2 1 1]);
colorbar;
colormap("hsv");
title('V1 orientation map');
%% ------------------------------------------------------------------------
%  3) Build demo Gabor filters for V1 and save
%     (Use the left/right orientation maps as inputs to the filter generator)
%  ------------------------------------------------------------------------
disp('Generating Gabor filters for V1...');

% Pull V1 RF/rate settings safely from struct p
V1_sigma   = p.V1.rf.sigma;         
V1_lambda  = p.V1.rf.lambda;
V1_rX_mean = p.V1.rates.mean_when_image;

[F, F_left, F_right] = physiological_properties.generate_Gabor_filter( ...
    thetaMap_V1_left, thetaMap_V1_right, ...
    'Sigma',  V1_sigma, ...
    'Lambda', V1_lambda, ...
    'rX',     V1_rX_mean ...
);

save(fullfile(dataDir, 'gaborFilters_V1.mat'), ...
     'F', 'F_left', 'F_right', 'thetaMap_V1_left', 'thetaMap_V1_right');

%% ------------------------------------------------------------------------
%  4) Generate demo weights (feedforward & recurrent) and save
%  ------------------------------------------------------------------------
disp('Generating demo connectivity weights...');

Ne = p.Ne; Ni = p.Ni; Nx = p.Nx;

% Feedforward and recurrent spatial scales (vectorized forms)
sigmaRX = p.conn.sigma.feedforward * ones(2,1); % [sigma_e; sigma_i]
sigmaRR = [ ...
    p.conn.sigma.recurrent_exc * ones(2,1), ... % E→{E,I}
    p.conn.sigma.recurrent_inh * ones(2,1)  ... % I→{E,I}
];

% Connection probabilities
Prr = [p.conn.pbar.ee, p.conn.pbar.ei; ...
       p.conn.pbar.ie, p.conn.pbar.ii];
Prx = [p.conn.pbar.eF; p.conn.pbar.iF];

% Tuning parameters
P_ts  = p.conn.excit.tuned_fraction * ones(2,1);
TS_th = p.conn.excit.tuned_cos_thresh;

[Wrr, Wrf] = physiological_properties.generate_weights_default( ...
    Ne, Ni, Nx, sigmaRX, sigmaRR, Prr, Prx, P_ts, TS_th, ...
    'ThetaV1',thetaMap_V1', 'ThetaE', thetaMap_V4MT');

save(fullfile(dataDir, 'demo_weight.mat'), ...
    'Wrr', 'Wrf', 'Ne', 'Ni', 'Nx', 'sigmaRX', 'sigmaRR', 'Prr', 'Prx', 'P_ts', 'TS_th');

%% ------------------------------------------------------------------------
%  5) Quick V1 firing-rate preview and plot
%  ------------------------------------------------------------------------
disp('Generating a quick V1 rate preview...');

out = utils_simulation.gen_V1_FiringRates(theta1/180*pi, theta2/180*pi, c1, c2);
if isfield(out, 'fr_left')
    fr = out.fr_left;   % pick one view for a simple preview
    save(fullfile(dataDir, 'demo_V1rate.mat'), 'fr');

    figure('Color','w','Name','V1 Firing Rate Preview');
    imagesc(reshape(fr, Nx_V1, Nx_V1/2)');  % 100 x 50 preview
    axis image; colorbar;
    title('Demo V1 rate (kHz)');
    pbaspect([1 2 1]);
else
    warning('gen_V1_FiringRates output missing field "fr_left". Skipping preview plot.');
end

%% ------------------------------------------------------------------------
%  6) Short simulation
%  ------------------------------------------------------------------------

% Prepare parameter overrides for the demo
weightfile = 'weightSpatRecTrans_sigmaRX_0d10_sigmaRR_0d20_Pts_0d15_0d15_tuning_th_0d60.mat';
param_changes = { ...
    'filename',  fullfile(dataDir, 'demo_simulation.mat'); ...
    'T',         T; ...
    'fr_fname',  fullfile(dataDir, 'demo_V1rate.mat'); ...
    'W_fname',   fullfile(dataDir, weightfile); ...
};

try
    spiking_simulation.spiking_sim(simulation_type, [], param_changes);
    disp('Simulation finished. Loading results...');
catch ME
    warning('Simulation failed: %s', ME.message);
end


%% ------------------------------------------------------------------------
%  7) Visualize MT population rate and a short raster (if available)
%  ------------------------------------------------------------------------

S = load(fullfile(dataDir, 'demo_simulation.mat'));  % expects s2 + param

s2    = S.s2;
param = S.params;

Ne_local = param.Ne; % fallback shape
MTfr = utils_analysis.spktime2count(s2, 1:Ne_local, param.T, 1, 1) / param.T * 1e3;

figure('Color','w','Name','MT Population Rate')

imagesc(reshape(MTfr, 200, 100)');
axis image; colorbar;
title('Demo MT population rate');
pbaspect([2 1 1]);


try
    figure('Color','w','Name','Raster Animation (short clip)');
    utils_analysis.raster2D_ani(s2, 1000, 2000, 100); % 1–2 s window, 100x100 layout
catch ME
    warning('Raster animation failed: %s', ME.message);
end
        

disp('All done! Explore saved outputs in the data/ and result/ folders.');

%% ------------------------------------------------------------------------
%  Local helpers
%  ------------------------------------------------------------------------
function ensureDir(p)
    if ~exist(p, 'dir'); mkdir(p); end
end

function assertAllExist(fqns)
    for i = 1:numel(fqns)
        if exist(fqns{i}, 'file') ~= 2
            error(['Required function not found on path: ', fqns{i}, ...
                   '\nCheck that your src/ folders are added and MEX files are compiled where needed.']);
        end
    end
end


function f = pickFunction(candidates)
%PICKFUNCTION  Return a function handle to the first existing candidate.
%   candidates: cellstr of names like {'spktime2count','utils_analysis.spktime2count'}
    f = [];
    for i = 1:numel(candidates)
        if exist(candidates{i}, 'file') == 2
            f = str2func(candidates{i});
            return;
        end
    end
end

