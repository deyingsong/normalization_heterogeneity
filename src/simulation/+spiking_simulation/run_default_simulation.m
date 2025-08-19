function run_default_simulation(params, options)
% RUN_DEFAULT_SIMULATION - Execute standard simulation
%
% H1 Line: Run default spiking network simulation
%
% Syntax:
%   run_default_simulation(params, options)
%
% Description:
%   Executes the standard simulation without additional noise or
%   special connectivity patterns.
%
% Inputs:
%   params  - Network parameters
%   options - Simulation options
%
% Example:
%   run_default_simulation(params, options);

    fprintf('Running default simulation...\n');
    
    % Load filters
    filters = utils_simulation.load_filters('V1filterRecSig0d2Lam0d6.mat');
    
    % Generate input spikes
    rng_seed = utils_simulation.get_rng_seed(options);
    load(params.fr_fname, 'fr');
    s1 = utils_simulation.gen_input_spikes(fr, params.T, filters.F, 450, params.sigma_n, params.tau_n);
    
    % Load or generate connectivity
    [Wrr, Wrf] = utils_simulation.gen_connectivity(params, options, rng_seed);
    
    % Initialize membrane potentials
    V0 = utils_simulation.initialize_membrane_potentials(params);
    params.V0 = V0;
    
    % Select recording neurons
    params.Irecord = utils_simulation.select_recording_neurons(params);
    
    % Run simulation
    tic;
    [s2, Isyn, Vm] = spiking_simulation.EIF_normalization_Default(s1, Wrf, Wrr, params);
    elapsed_time = toc;
    
    % Analyze results
    results = utils_simulation.analyze_simulation_results(s2, params, elapsed_time);
    
    % Save results
    if options.save
        utils_simulation.save_simulation_results(params.filename, s2, s1, results, params, options);
    end
end
