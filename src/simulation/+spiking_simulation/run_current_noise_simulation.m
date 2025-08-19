function run_current_noise_simulation(params, options)
% RUN_CURRENT_NOISE_SIMULATION - Run with current noise injection
%
% H1 Line: Execute simulation with additive current noise
%
% Syntax:
%   run_current_noise_simulation(params, options)
%
% Description:
%   Runs simulation with white noise current injection to each neuron
%   to study noise effects on network dynamics.
%
% Inputs:
%   params  - Network parameters
%   options - Simulation options
%
% Example:
%   run_current_noise_simulation(params, options);

    fprintf('Running simulation with current noise...\n');
    fprintf('Noise level: sigma = %.2f\n', params.sigma_current);
    
    % Standard initialization
    filters = load_filters('V1filterRecSig0d2Lam0d6.mat');
    rng_seed = get_rng_seed(options);
    s1 = generate_input_spikes(params, options, filters, rng_seed);
    [Wrr, Wrf] = generate_connectivity(params, options, rng_seed);
    
    % Set random seed for noise generation
    rng(rng_seed);
    
    % Initialize and run
    V0 = initialize_membrane_potentials(params);
    params.V0 = V0;
    params.Irecord = select_recording_neurons(params);
    
    tic;
    [s2, Isyn, Vm] = EIF1DRFfastslowSynAtttSpatRec_CurrentNoise(...
        s1, Wrf, Wrr, params);
    elapsed_time = toc;
    
    % Analyze and save
    results = analyze_simulation_results(s2, params, elapsed_time);
    
    if options.save
        save_simulation_results(options.filename, s2, results, params, options);
        if ~isempty(Isyn) || ~isempty(Vm)
            save(options.filename, 'Isyn', 'Vm', '-append');
        end
    end
end
