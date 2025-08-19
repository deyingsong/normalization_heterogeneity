function out = generate_physiological_properties()
%GENERATE_PHYSIOLOGICAL_PROPERTIES Build maps, filters, and weights, then save to /data.
% Syntax
%   out = generate_physiological_properties()
%
% Description
%   - Adds src/simulation to path and loads network/constants via neuronConstants()
%   - Generates orientation maps for V1 (Nx=100, split into two 50x50 fields)
%     and for V4/MT (Nx=200, split into two 100x100 fields)
%   - Builds Gabor filters from the two V1 subfields
%   - Generates default recurrent (Wrr) and feedforward (Wrf) weights
%   - Saves named .mat files under <repoRoot>/data and returns their paths
%
% Output (struct)
%   out.paths.thetaMaps.V1.all, out.paths.thetaMaps.V1.left, out.paths.thetaMaps.V1.right
%   out.paths.thetaMaps.V4MT.all, out.paths.thetaMaps.V4MT.left, out.paths.thetaMaps.V4MT.right
%   out.paths.gabor.V1
%   out.paths.weights.default
%
% Notes
%   - Requires: att_ori_map.m, physiological_properties.generate_Gabor_filter.m,
%               generate_weights_default.m, neuronConstants.m
%   - All variables are saved with descriptive names.

  %% Resolve repo paths and add sim sources
  thisFile  = mfilename('fullpath');
  repoRoot  = fileparts(fileparts(thisFile));   % up from sim_script → repo root
  simSrcDir = fullfile(repoRoot, 'src', 'simulation');
  if ~exist(simSrcDir, 'dir')
      error('Simulation source folder not found: %s', simSrcDir);
  end
  addpath(simSrcDir);

  dataDir = fullfile(repoRoot, 'data');
  if ~exist(dataDir, 'dir'), mkdir(dataDir); end

  %% Load constants / params
  p = neuronConstants();  % 'p' used below

  %% -------------------------
  %  V1 orientation map (Nx=100) and subfields (50×50 each)
  %  -------------------------
  Nx_V1 = 100;
  Kc    = 10 * 2*pi;  % as in your snippet

  thetaMap_V1 = physiological_properties.generate_ori_map(Nx_V1, Kc);   % size 100×100
  thetaMap_V1_topHalf  = thetaMap_V1(1:50, :);
  thetaMap_V1_left     = thetaMap_V1_topHalf(:, 1:50);
  thetaMap_V1_right    = thetaMap_V1_topHalf(:, 51:100);

  save(fullfile(dataDir, 'thetaMap_V1_all.mat'),       'thetaMap_V1');
  save(fullfile(dataDir, 'thetaMap_V1_left.mat'),      'thetaMap_V1_left');
  save(fullfile(dataDir, 'thetaMap_V1_right.mat'),     'thetaMap_V1_right');

  %% -------------------------
  %  V4/MT orientation map (Nx=200) and subfields (100×100 each)
  %  -------------------------
  Nx_V4MT = 200;
  thetaMap_V4MT = physiological_properties.generate_ori_map(Nx_V4MT, Kc);   % size 200×200

  % Split into two 100×100 fields to mirror the V1 handling
  thetaMap_V4MT_topHalf = thetaMap_V4MT(1:100, :);
  thetaMap_V4MT_left    = thetaMap_V4MT_topHalf(:, 1:100);
  thetaMap_V4MT_right   = thetaMap_V4MT_topHalf(:, 101:200);

  save(fullfile(dataDir, 'thetaMap_V4MT_all.mat'),     'thetaMap_V4MT');
  save(fullfile(dataDir, 'thetaMap_V4MT_left.mat'),    'thetaMap_V4MT_left');
  save(fullfile(dataDir, 'thetaMap_V4MT_right.mat'),   'thetaMap_V4MT_right');

  %% -------------------------
  %  Gabor filters for the two V1 subfields
  %  -------------------------
  % Use the two 50×50 V1 fields as orientation maps for filter construction:
  ori_map_V11 = thetaMap_V1_left;
  ori_map_V12 = thetaMap_V1_right;

  % Allow fields in 'p' to drive parameters (names taken from your snippet)
  % p.V1.rf.sigma, p.V1.rf.lambda, p.V1.rates.mean_when_image
  [F, F_left, F_right] = physiological_properties.generate_Gabor_filter( ...
      ori_map_V11, ori_map_V12, ...
      'Sigma',  getfield(p, 'V1').rf.sigma, ...        
      'Lambda', getfield(p, 'V1').rf.lambda, ...       
      'rX',     getfield(p, 'V1').rates.mean_when_image ... 
  );

  save(fullfile(dataDir, 'gaborFilters_V1.mat'), ...
       'F', 'F_left', 'F_right', 'ori_map_V11', 'ori_map_V12');

  %% -------------------------
  %  Default connectivity weights
  %  -------------------------
  Ne = p.Ne;
  Ni = p.Ni;
  Nx = p.Nx;

  % Feedforward and recurrent spatial scales
  sigmaRX = p.conn.sigma.feedforward * ones(2,1); % [sigma_e; sigma_i] for X→{E,I}
  sigmaRR = [ ...
      p.conn.sigma.recurrent_exc * ones(2,1), ... % E→{E,I}
      p.conn.sigma.recurrent_inh * ones(2,1) ...  % I→{E,I}
  ];

  % Connection probabilities (rows: to E; to I) × (from E, from I)
  % and feedforward to {E,I}
  Prr = [p.conn.pbar.ee, p.conn.pbar.ei; ...
         p.conn.pbar.ie, p.conn.pbar.ii];
  Prx = [p.conn.pbar.eF; p.conn.pbar.iF];

  % Tuning settings (fraction of tuned E synapses and cosine threshold)
  P_ts  = p.conn.excit.tuned_fraction*ones(2,1);
  TS_th = p.conn.excit.tuned_cos_thresh;  % fixed name and source

  % Generate weights
  [Wrr, Wrf] = physiological_properties.generate_weights_default( ...
      Ne, Ni, Nx, sigmaRX, sigmaRR, Prr, Prx, P_ts, TS_th);

  save(fullfile(dataDir, 'weights_default.mat'), ...
       'Wrr', 'Wrf', 'Ne', 'Ni', 'Nx', 'sigmaRX', 'sigmaRR', 'Prr', 'Prx', 'P_ts', 'TS_th');

  %% Output paths for convenience
  out = struct();
  out.paths.thetaMaps.V1.all   = fullfile(dataDir, 'thetaMap_V1_all.mat');
  out.paths.thetaMaps.V1.left  = fullfile(dataDir, 'thetaMap_V1_left.mat');
  out.paths.thetaMaps.V1.right = fullfile(dataDir, 'thetaMap_V1_right.mat');

  out.paths.thetaMaps.V4MT.all   = fullfile(dataDir, 'thetaMap_V4MT_all.mat');
  out.paths.thetaMaps.V4MT.left  = fullfile(dataDir, 'thetaMap_V4MT_left.mat');
  out.paths.thetaMaps.V4MT.right = fullfile(dataDir, 'thetaMap_V4MT_right.mat');

  out.paths.gabor.V1    = fullfile(dataDir, 'gaborFilters_V1.mat');
  out.paths.weights.default = fullfile(dataDir, 'weights_default.mat');

  fprintf('Saved:\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n', ...
      out.paths.thetaMaps.V1.all, out.paths.thetaMaps.V1.left, out.paths.thetaMaps.V1.right, ...
      out.paths.thetaMaps.V4MT.all, out.paths.thetaMaps.V4MT.left, out.paths.thetaMaps.V4MT.right, ...
      out.paths.gabor.V1);
  fprintf('  %s\n', out.paths.weights.default);
end
