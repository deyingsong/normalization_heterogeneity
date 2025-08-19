function run_broad_weight_simulation(params, options)
% RUN_BROAD_WEIGHT_SIMULATION - Run simulation with broad weight distribution
%
% H1 Line: Execute simulation with gamma-distributed synaptic weights
%
% Syntax:
%   run_broad_weight_simulation(params, options)
%
% Description:
%   Runs simulation where synaptic weights follow gamma distributions
%   for more realistic weight heterogeneity.
%
% Inputs:
%   params  - Network parameters
%   options - Simulation options
%
% Example:
%   run_broad_weight_simulation(params, options);

    fprintf('Running broad weight distribution simulation...\n');
    
    % Check for required parameters
    required_params = {'bee', 'bei', 'bie', 'bii', 'bex', 'bix'};
    for i = 1:length(required_params)
        if ~isfield(params, required_params{i})
            error('SpikingSim:MissingParameter', ...
                  'Broad weight simulation requires parameter: %s', ...
                  required_params{i});
        end
    end
    
    % Load filters
    filters = load_filters('V1filterRecSig0d2Lam0d6.mat');
    
    % Generate input spikes
    rng_seed = get_rng_seed(options);
    s1 = generate_input_spikes(params, options, filters, rng_seed);
    
    % Load connectivity
    [Wrr, Wrf] = generate_connectivity(params, options, rng_seed);
    
    % Generate weight matrices with gamma distribution
    [Jrr, Jrf] = generate_gamma_weights(params, Wrr, Wrf);
    
    % Initialize and run simulation
    V0 = initialize_membrane_potentials(params);
    params.V0 = V0;
    params.Irecord = select_recording_neurons(params);
    
    tic;
    [s2, Isyn, Vm] = EIF1DRFfastslowSynAtttSpatRec(s1, Wrf, Wrr, Jrf, Jrr, params);
    elapsed_time = toc;
    
    % Analyze and save results
    results = analyze_simulation_results(s2, params, elapsed_time);
    
    if options.save
        save_simulation_results(options.filename, s2, results, params, options);
    end
end

function [Jrr, Jrf] = generate_gamma_weights(params, Wrr, Wrf)
% GENERATE_GAMMA_WEIGHTS - Create gamma-distributed synaptic weights
%
% H1 Line: Generate synaptic weight matrices with gamma distributions
%
% Syntax:
%   [Jrr, Jrf] = generate_gamma_weights(params, Wrr, Wrf)
%
% Description:
%   Creates synaptic weight matrices where weights follow gamma
%   distributions specified by shape parameters.
%
% Inputs:
%   params - Network parameters including shape parameters
%   Wrr    - Recurrent connectivity
%   Wrf    - Feedforward connectivity
%
% Outputs:
%   Jrr - Recurrent weight matrix
%   Jrf - Feedforward weight matrix
%
% Example:
%   [Jrr, Jrf] = generate_gamma_weights(params, Wrr, Wrf);

    % Initialize weight matrices
    Jrf = zeros(size(Wrf));
    Jrr = zeros(size(Wrr));
    
    % Generate feedforward weights
    Jex_samples = gamrnd(params.Jx(1)/params.bex, params.bex, ...
                        [params.Kx(1)*params.Nx, 1]);
    Jix_samples = gamrnd(params.Jx(1)/params.bix, params.bix, ...
                        [params.Kx(2)*params.Nx, 1]);
    
    % Fill feedforward weight matrix
    [Jrf] = fill_weight_matrix(Jrf, Jex_samples, Jix_samples, ...
                              params.Kx, params.Nx);
    
    % Generate recurrent weights
    Jee_samples = gamrnd(params.Jr(1,1)/params.bee, params.bee, ...
                        [params.Kr(1,1)*params.Ne, 1]);
    Jie_samples = gamrnd(params.Jr(2,1)/params.bie, params.bie, ...
                        [params.Kr(2,1)*params.Ne, 1]);
    Jei_samples = gamrnd(-params.Jr(1,2)/params.bei, params.bei, ...
                        [params.Kr(1,2)*params.Ni, 1]);
    Jii_samples = gamrnd(-params.Jr(2,2)/params.bii, params.bii, ...
                        [params.Kr(2,2)*params.Ni, 1]);
    
    % Fill recurrent weight matrix
    [Jrr] = fill_recurrent_weights(Jrr, Jee_samples, Jie_samples, ...
                                   Jei_samples, Jii_samples, params);
end

function matrix = fill_weight_matrix(matrix, exc_samples, inh_samples, K, N)
% FILL_WEIGHT_MATRIX - Fill weight matrix with samples
%
% H1 Line: Populate weight matrix with excitatory and inhibitory samples
%
% Syntax:
%   matrix = fill_weight_matrix(matrix, exc_samples, inh_samples, K, N)
%
% Description:
%   Fills a weight matrix with excitatory and inhibitory weight samples
%   in the appropriate positions.
%
% Inputs:
%   matrix      - Empty weight matrix to fill
%   exc_samples - Excitatory weight samples
%   inh_samples - Inhibitory weight samples
%   K          - Out-degree vector [K_exc; K_inh]
%   N          - Number of presynaptic neurons
%
% Outputs:
%   matrix - Filled weight matrix
%
% Example:
%   matrix = fill_weight_matrix(zeros(1000,1), exc_w, inh_w, [10;5], 50);

    m = K(1); % Excitatory out-degree
    n = K(2); % Inhibitory out-degree
    
    % Calculate indices for excitatory connections
    exc_idx = reshape(bsxfun(@plus, (0:N-1)' * (m+n), 1:m), [], 1);
    
    % Calculate indices for inhibitory connections
    inh_idx = reshape(bsxfun(@plus, (0:N-1)' * (m+n), m+1:m+n), [], 1);
    
    % Fill matrix
    matrix(exc_idx) = exc_samples;
    matrix(inh_idx) = inh_samples;
end


function Jrr = fill_recurrent_weights(Jrr, Jee, Jie, Jei, Jii, params)
% FILL_RECURRENT_WEIGHTS - Fill recurrent weight matrix
%
% H1 Line: Populate recurrent connectivity with weight samples
%
% Syntax:
%   Jrr = fill_recurrent_weights(Jrr, Jee, Jie, Jei, Jii, params)
%
% Description:
%   Fills recurrent weight matrix with E-E, I-E, E-I, and I-I weights
%   in the appropriate block structure.
%
% Inputs:
%   Jrr    - Empty recurrent weight matrix
%   Jee    - E->E weight samples
%   Jie    - I->E weight samples
%   Jei    - E->I weight samples
%   Jii    - I->I weight samples
%   params - Network parameters
%
% Outputs:
%   Jrr - Filled recurrent weight matrix
%
% Example:
%   Jrr = fill_recurrent_weights(Jrr, Jee, Jie, Jei, Jii, params);

    % Fill E->E and I->E connections
    offset = 0;
    Jrr = fill_weight_block(Jrr, Jee, Jie, params.Kr(:,1), ...
                           params.Ne, offset);
    
    % Fill E->I and I->I connections
    offset = params.Ne * sum(params.Kr(:,1));
    Jrr = fill_weight_block(Jrr, Jei, Jii, params.Kr(:,2), ...
                           params.Ni, offset);
end


function matrix = fill_weight_block(matrix, samples1, samples2, K, N, offset)
% FILL_WEIGHT_BLOCK - Helper to fill weight matrix blocks
%
% H1 Line: Fill a block of the weight matrix
%
% Syntax:
%   matrix = fill_weight_block(matrix, samples1, samples2, K, N, offset)
%
% Inputs:
%   matrix   - Weight matrix to modify
%   samples1 - First set of weight samples
%   samples2 - Second set of weight samples
%   K       - Out-degrees [K1; K2]
%   N       - Number of neurons
%   offset  - Starting index offset
%
% Outputs:
%   matrix - Updated weight matrix

    m = K(1);
    n = K(2);
    
    idx1 = reshape(bsxfun(@plus, (0:N-1)' * (m+n), 1:m), [], 1) + offset;
    idx2 = reshape(bsxfun(@plus, (0:N-1)' * (m+n), m+1:m+n), [], 1) + offset;
    
    matrix(idx1) = samples1;
    matrix(idx2) = samples2;
end
