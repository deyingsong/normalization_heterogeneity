function s1 = generate_input_spikes(params, options, filters, rng_seed)
% GENERATE_INPUT_SPIKES - Generate input spike trains
%
% H1 Line: Create input spike trains for feedforward layer
%
% Syntax:
%   s1 = generate_input_spikes(params, options, filters, rng_seed)
%
% Description:
%   Generates spike trains for the feedforward input layer based on
%   firing rates and specified noise parameters.
%
% Inputs:
%   params   - Network parameters
%   options  - Simulation options
%   filters  - Filter matrix for input generation
%   rng_seed - Random seed for reproducibility
%
% Outputs:
%   s1 - Spike train matrix [time, neuron_id]
%
% Example:
%   s1 = generate_input_spikes(params, options, F, 12345);

    % Set random seed
    if nargin >= 4 && ~isempty(rng_seed)
        rng(rng_seed);
    end
    
    % Load firing rates
    if options.loadfr1
        if ~exist(options.fr_fname1, 'file')
            error('SpikingSim:FileNotFound', ...
                  'Firing rate file not found: %s', options.fr_fname1);
        end
        
        fr_data = load(options.fr_fname1, 'fr');
        fr = fr_data.fr;
        
        % Validate firing rates
        if any(fr < 0)
            error('SpikingSim:InvalidFiringRates', ...
                  'Negative firing rates detected');
        end
        
        % Generate spikes with noise
        s1 = attt_genXspk_noise(fr, params.T, filters, params.NI, ...
                               params.sigma_n, params.tau_n);
    else
        % Generate uniform random spikes if no rates provided
        warning('SpikingSim:NoFiringRates', ...
                'No firing rates provided, using random spikes');
        s1 = generate_random_spikes(params);
    end
    
    % Validate output
    if isempty(s1)
        error('SpikingSim:EmptySpikes', ...
              'No input spikes generated');
    end
end

function s1 = generate_random_spikes(params)
% GENERATE_RANDOM_SPIKES - Create random spike trains
%
% H1 Line: Generate random Poisson spike trains
%
% Syntax:
%   s1 = generate_random_spikes(params)
%
% Description:
%   Creates random spike trains with Poisson statistics for testing
%   or when no firing rate data is available.
%
% Inputs:
%   params - Network parameters
%
% Outputs:
%   s1 - Random spike train matrix
%
% Example:
%   s1 = generate_random_spikes(params);

    base_rate = 10; % Hz
    n_spikes = poissrnd(base_rate * params.T / 1000 * params.Nx);
    
    s1 = zeros(2, n_spikes);
    s1(1, :) = sort(rand(1, n_spikes) * params.T);
    s1(2, :) = randi(params.Nx, 1, n_spikes);
end
