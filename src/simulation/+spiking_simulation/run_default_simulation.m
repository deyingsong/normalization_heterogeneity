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
    filters = load_filters('V1filterRecSig0d2Lam0d6.mat');
    
    % Generate input spikes
    rng_seed = get_rng_seed(options);
    s1 = generate_input_spikes(params, options, filters, rng_seed);
    
    % Load or generate connectivity
    [Wrr, Wrf] = generate_connectivity(params, options, rng_seed);
    
    % Initialize membrane potentials
    V0 = initialize_membrane_potentials(params);
    params.V0 = V0;
    
    % Select recording neurons
    params.Irecord = select_recording_neurons(params);
    
    % Run simulation
    tic;
    [s2, Isyn, Vm] = EIF1DRFfastslowSynAtttSpatRec(s1, Wrf, Wrr, params);
    elapsed_time = toc;
    
    % Analyze results
    results = analyze_simulation_results(s2, params, elapsed_time);
    
    % Save results
    if options.save
        save_simulation_results(options.filename, s2, results, params, options);
    end
end
