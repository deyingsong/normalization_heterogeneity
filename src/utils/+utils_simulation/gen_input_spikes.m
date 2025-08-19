function sx = generate_input_spikes(fr, T, varargin)
% GENERATE_INPUT_SPIKES - Generate spike trains from firing rates
%
% H1 Line: Create Poisson spike trains with optional noise modulation
%
% Syntax:
%   sx = generate_input_spikes(fr, T)
%   sx = generate_input_spikes(fr, T, F, NI, sigma_n, tau_n)
%   sx = generate_input_spikes(fr, T, 'NoiseParams', params)
%
% Description:
%   Generates spike trains for input neurons based on firing rates.
%   Can include noise modulation through filter matrix F. Uses Poisson
%   process for spike generation with time-varying rates.
%
% Inputs:
%   fr      - Vector of baseline firing rates (Hz) for each neuron
%   T       - Total simulation time (ms)
%   F       - Filter matrix (optional, for noise modulation)
%   NI      - Number of noise inputs (optional)
%   sigma_n - Noise standard deviation (optional)
%   tau_n   - Noise correlation time (optional)
%
% Outputs:
%   sx - Spike train matrix (2 x nspikes)
%        sx(1,:) - spike times (ms)
%        sx(2,:) - neuron indices
%
% Example:
%   % Simple Poisson spikes
%   sx = generate_input_spikes(ones(100,1)*10, 1000);
%   
%   % With noise modulation
%   sx = generate_input_spikes(fr, 1000, F, 450, 3.5, 40);
%
% See also: generate_fisher_info_spikes

    % Parse inputs
    [F, NI, sigma_n, tau_n, has_noise] = parse_spike_inputs(varargin{:});
    
    % Validate inputs
    validate_spike_parameters(fr, T);
    
    % Generate spikes
    if has_noise
        sx = generate_noisy_spikes(fr, T, F, NI, sigma_n, tau_n);
    else
        sx = generate_poisson_spikes(fr, T);
    end
end

function [F, NI, sigma_n, tau_n, has_noise] = parse_spike_inputs(varargin)
% PARSE_SPIKE_INPUTS - Parse optional inputs for spike generation
%
% H1 Line: Parse variable inputs for spike generation function

    has_noise = false;
    F = [];
    NI = [];
    sigma_n = [];
    tau_n = [];
    
    if nargin == 4
        F = varargin{1};
        NI = varargin{2};
        sigma_n = varargin{3};
        tau_n = varargin{4};
        has_noise = true;
    elseif nargin == 2 && isstruct(varargin{2})
        % Support for parameter struct
        params = varargin{2};
        if isfield(params, 'F')
            F = params.F;
            NI = params.NI;
            sigma_n = params.sigma_n;
            tau_n = params.tau_n;
            has_noise = true;
        end
    elseif nargin > 0
        error('Spikes:InvalidInputs', ...
              'Invalid number of input arguments');
    end
end

function validate_spike_parameters(fr, T)
% VALIDATE_SPIKE_PARAMETERS - Check spike generation parameters
%
% H1 Line: Validate firing rates and time duration

    % Check firing rates
    if ~isvector(fr) || ~isnumeric(fr)
        error('Spikes:InvalidFR', 'Firing rates must be a numeric vector');
    end
    
    if any(fr < 0)
        error('Spikes:NegativeFR', 'Firing rates cannot be negative');
    end
    
    if any(fr > 1000)
        warning('Spikes:HighFR', ...
                'Very high firing rates detected (max: %.1f Hz)', max(fr));
    end
    
    % Check time
    if ~isscalar(T) || T <= 0
        error('Spikes:InvalidT', 'T must be a positive scalar');
    end
    
    % Ensure column vector
    if isrow(fr)
        fr = fr';
    end
end

function sx = generate_poisson_spikes(fr, T)
% GENERATE_POISSON_SPIKES - Generate simple Poisson spike trains
%
% H1 Line: Create Poisson spikes with constant rates
%
% Syntax:
%   sx = generate_poisson_spikes(fr, T)
%
% Description:
%   Generates spike trains using a Poisson process with constant
%   firing rates for each neuron.
%
% Inputs:
%   fr - Firing rate vector (Hz)
%   T  - Duration (ms)
%
% Outputs:
%   sx - Spike train matrix

    Nx = length(fr);
    dt = 1;  % 1 ms time step
    
    % Pre-allocate with expected size
    expected_spikes = sum(fr) * T / 1000;
    sx = zeros(2, round(expected_spikes * 1.5));  % 50% safety margin
    
    nspks = 0;
    
    % Generate spikes for each time step
    for t = dt:dt:T
        % Poisson spike generation
        spk = rand(Nx, 1) < fr * dt / 1000;  % Convert Hz to probability
        
        ns_temp = nnz(spk);
        if ns_temp > 0
            % Store spike times and neuron indices
            sx(1, nspks+(1:ns_temp)) = t;
            sx(2, nspks+(1:ns_temp)) = find(spk);
            nspks = nspks + ns_temp;
        end
    end
    
    % Trim to actual size and sort
    sx = sx(:, 1:nspks);
    [~, J] = sort(sx(1, :));
    sx = sx(:, J);
end

function sx = generate_noisy_spikes(fr, T, F, NI, sigma_n, tau_n)
% GENERATE_NOISY_SPIKES - Generate spikes with noise modulation
%
% H1 Line: Create spike trains with time-varying noisy rates
%
% Syntax:
%   sx = generate_noisy_spikes(fr, T, F, NI, sigma_n, tau_n)
%
% Description:
%   Generates spike trains where firing rates are modulated by
%   filtered noise inputs, creating more realistic variability.
%
% Inputs:
%   fr      - Baseline firing rates (Hz)
%   T       - Duration (ms)
%   F       - Filter matrix
%   NI      - Number of noise inputs
%   sigma_n - Noise standard deviation
%   tau_n   - Noise time constant
%
% Outputs:
%   sx - Spike train matrix

    Nx = length(fr);
    dt = 1;  % 1 ms time step
    
    % Validate filter matrix dimensions
    if size(F, 1) ~= Nx || size(F, 2) ~= NI
        error('Spikes:FilterDimension', ...
              'Filter matrix size (%dx%d) incompatible with Nx=%d, NI=%d', ...
              size(F, 1), size(F, 2), Nx, NI);
    end
    
    % Pre-allocate
    expected_spikes = sum(fr) * T / 1000 * 1.2;  % Extra for noise
    sx = zeros(2, round(expected_spikes * 1.5));
    nspks = 0;
    
    % Initialize noise with stationary distribution
    noise = sigma_n / sqrt(2 * tau_n) * randn(NI, 1);
    
    % Generate spikes with time-varying rates
    for t = dt:dt:T
        % Update noise (Ornstein-Uhlenbeck process)
        noise = noise + (-noise * dt + sigma_n * randn(NI, 1) * sqrt(dt)) / tau_n;
        
        % Calculate instantaneous firing rates
        FR = fr + F * noise;
        
        % Rectify negative rates
        FR(FR < 0) = 0;
        
        % Generate spikes
        spk = rand(Nx, 1) < FR * dt / 1000;
        
        ns_temp = nnz(spk);
        if ns_temp > 0
            sx(1, nspks+(1:ns_temp)) = t;
            sx(2, nspks+(1:ns_temp)) = find(spk);
            nspks = nspks + ns_temp;
        end
    end
    
    % Trim and sort
    sx = sx(:, 1:nspks);
    [~, J] = sort(sx(1, :));
    sx = sx(:, J);
end
