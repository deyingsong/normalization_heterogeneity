function run_gaussian_noise_simulation(params, options)
% RUN_GAUSSIAN_NOISE_SIMULATION - Run with Gaussian noise modulation
%
% H1 Line: Execute simulation with Gaussian noise on inputs
%
% Syntax:
%   run_gaussian_noise_simulation(params, options)
%
% Description:
%   Runs simulation with temporally correlated Gaussian noise
%   modulating the input currents.
%
% Inputs:
%   params  - Network parameters
%   options - Simulation options
%
% Example:
%   run_gaussian_noise_simulation(params, options);

    fprintf('Running simulation with Gaussian noise modulation...\n');
    
    % Check required parameters
    if ~isfield(params, 'tau_L') || ~isfield(params, 'sigma_s')
        error('SpikingSim:MissingNoiseParams', ...
              'Gaussian noise requires tau_L and sigma_s parameters');
    end
    
    % Standard initialization
    filters = utils_simulation.load_filters('V1filterRecSig0d2Lam0d6.mat');
    rng_seed = utils_simulation.get_rng_seed(options);
    s1 = utils_simulation.gen_input_spikes(params, options, filters, rng_seed);
    [Wrr, Wrf] = utils_simulation.gen_connectivity(params, options, rng_seed);
    
    % Generate Gaussian noise signal
    L = utils_simulation.gen_Gaussian_noise(params.T, 1, params.tau_L, params.sigma_s);
    if ~iscolumn(L)
        L = L';
    end
    
    % Validate noise signal
    if length(L) ~= params.T/params.dt
        warning('SpikingSim:NoiseLengthMismatch', ...
                'Noise signal length does not match simulation duration');
    end
    
    % Initialize and run
    V0 = utils_simulation.initialize_membrane_potentials(params);
    params.V0 = V0;
    params.Irecord = utils_simulation.select_recording_neurons(params);
    
    tic;
    [s2, Isyn, Vm] = spiking_simulation.EIF_normalization_CurrentGaussianNoise(...
        s1, Wrf, Wrr, params, L);
    elapsed_time = toc;
    
    % Analyze and save
    results = utils_simulation.analyze_simulation_results(s2, params, elapsed_time);
    
    if options.save
        utils_simulation.save_simulation_results(options.filename, s2, results, params, options);
        save(options.filename, 'L', '-append');
    end
end
