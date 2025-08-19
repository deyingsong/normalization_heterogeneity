function run_fisher_info_simulation(params, options, match_indegree)
% RUN_FISHER_INFO_SIMULATION - Run simulation for Fisher information analysis
%
% H1 Line: Execute simulation optimized for Fisher information calculation
%
% Syntax:
%   run_fisher_info_simulation(params, options, match_indegree)
%
% Description:
%   Runs simulation with ON/OFF cell separation and spike count
%   recording for Fisher information analysis.
%
% Inputs:
%   params         - Network parameters
%   options        - Simulation options
%   match_indegree - Boolean flag for in-degree matching
%
% Example:
%   run_fisher_info_simulation(params, options, false);

    fprintf('Running Fisher information simulation...\n');
    if match_indegree
        fprintf('Using in-degree matched connectivity\n');
    end
    
    thisFile = mfilename('fullpath');
    srcDir = fileparts(fileparts(thisFile));
    utilsFolder = fullfile(srcDir, 'utils', '+utils_simulation');
    utils_name = fullfile(srcDir, 'utils', '+utils_analysis', 'spktime2count.c');
    mex('-silent', '-outdir', utilsFolder, utils_name);

    % Check for contrast percentage parameter
    if ~isfield(params, 'ConPerc')
        error('SpikingSim:MissingConPerc', ...
              'Fisher info simulation requires ConPerc parameter');
    end
    
    % Load and prepare filters for ON/OFF cells
    filters = utils_simulation.load_filters('V1filterRecSig0d2Lam0d6.mat');
    [F1, F2] = prepare_filters(filters.F);
    
    % Generate input spikes with ON/OFF structure
    rng_seed = utils_simulation.get_rng_seed(options);
    fr_data = load(options.fr_fname1, 'fr');
    s1 = utils_simulation.gen_fisher_info_spikes(fr_data.fr, params.T, F1, F2, ...
                                     params.NI, params.sigma_n, ...
                                     params.tau_n, params.ConPerc);
    
    % Load connectivity
    if match_indegree
        [Wrr, Wrf, degree_params] = utils_simulation.load_indegree_connectivity(options, params);
        params = utils_simulation.merge_structs(params, degree_params);
    else
        [Wrr, Wrf] = utils_simulation.gen_connectivity(params, options, rng_seed);
    end
    
    % Initialize
    V0 = utils_simulation.initialize_membrane_potentials(params);
    params.V0 = V0;
    
    % Set more recording neurons for Fisher info
    params.nrecordE0 = 100;
    params.nrecordI0 = 100;
    params.Irecord = utils_simulation.select_recording_neurons(params);
    
    % Run simulation
    tic;
    if match_indegree
        [s2, ~, ~] = spiking_simulation.EIF_normalization_MatchInDegree(...
            s1, Wrf, Wrr, params);
    else
        [s2, ~, ~] = spiking_simulation.EIF_normalization_Default(...
            s1, Wrf, Wrr, params);
    end
    elapsed_time = toc;
    
    % Process for Fisher information
    s2 = s2(:, s2(2,:) ~= 0);
    MTE2 = utils_simulation.spktime2count(s2, 1:params.T, 100, 200, 1);
    
    % Process spike counts for analysis windows
    timeinds = reshape(((3:40)-1)*5+(1:2)', [], 1);
    MTE2_process = MTE2(:, timeinds);
    MTE2_process = MTE2_process(:, 1:2:75) + MTE2_process(:, 2:2:76);
    MTE2_process = int8(MTE2_process');
    
    % Calculate results
    results = utils_simulation.analyze_simulation_results(s2, params, elapsed_time);
    
    % Save
    if options.save
        save(options.filename, 's2', 'MTE2', 'MTE2_process', ...
             'results', 'params', '-v7.3');
    end
end

function [F1, F2] = prepare_filters(F)
% PREPARE_ON_OFF_FILTERS - Split filters for V1 area 1 and V1 area 2 channels
%
% H1 Line: Separate visual filters into ON and OFF pathways
%
% Syntax:
%   [F1, F2] = prepare_on_off_filters(F)
%
% Description:
%   Splits the filter matrix into V1 area 1 and area 2 cell populations
%   for modeling separate visual pathways.
%
% Inputs:
%   F - Original filter matrix (100x100 x channels)
%
% Outputs:
%   F1 - ON cell filters (50x50 x channels/2)
%   F2 - OFF cell filters (50x50 x channels/2)
%
% Example:
%   [F1, F2] = prepare_on_off_filters(F);

    % Validate input
    if size(F, 1) ~= 100 || size(F, 2) ~= 100
        error('SpikingSim:InvalidFilterSize', ...
              'Filter matrix must be 100x100xN');
    end
    
    % Initialize output matrices
    F1 = zeros(50*50, 225);
    F2 = zeros(50*50, 225);
    
    % Split into ON and OFF populations
    for i = 1:50
        for j = 1:50
            for k = 1:225
                F1((i-1)*50+j, k) = F((i-1)*100+j, k);
                F2((i-1)*50+j, k) = F((i-1)*100+j+50, k+225);
            end
        end
    end
end
