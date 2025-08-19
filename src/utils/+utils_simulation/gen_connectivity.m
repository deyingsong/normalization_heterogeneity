function [Wrr, Wrf] = gen_connectivity(params, options, rng_seed)
% GENERATE_CONNECTIVITY - Create network connectivity matrices
%
% H1 Line: Generate connection matrices for neural network
%
% Syntax:
%   [Wrr, Wrf] = generate_connectivity(params, options, rng_seed)
%
% Description:
%   Generates connectivity matrices for recurrent and feedforward
%   connections based on specified parameters and options.
%
% Inputs:
%   params   - Network parameters structure
%   options  - Simulation options structure
%   rng_seed - Random seed for reproducibility
%
% Outputs:
%   Wrr - Recurrent connectivity matrix
%   Wrf - Feedforward connectivity matrix
%
% Example:
%   [Wrr, Wrf] = generate_connectivity(params, options, 12345);

    % Set random seed for reproducibility
    if nargin >= 3 && ~isempty(rng_seed)
        rng(rng_seed, 'combRecursive');
    end
    
    if options.useWfile
        % Load existing connectivity
        if ~exist(params.W_fname, 'file')
            error('SpikingSim:FileNotFound', ...
                  'Weight file not found: %s', params.W_fname);
        end

        loaded_data = load(params.W_fname, 'Wrr', 'Wrf');
        Wrr = loaded_data.Wrr;
        Wrf = loaded_data.Wrf;

        fprintf('Loaded connectivity from: %s\n', params.W_fname);
    else
        % Generate new connectivity
        % This would call the actual connectivity generation function
        % Placeholder for actual implementation
        error('SpikingSim:NotImplemented', ...
              'Connectivity generation not yet implemented');
    end
    
    % Validate connectivity matrices
    validate_connectivity(Wrr, Wrf, params);
end

function validate_connectivity(Wrr, Wrf, params)
% VALIDATE_CONNECTIVITY - Check connectivity matrix validity
%
% H1 Line: Validate connectivity matrices for consistency
%
% Syntax:
%   validate_connectivity(Wrr, Wrf, params)
%
% Description:
%   Checks that connectivity matrices have correct dimensions and
%   properties for the specified network parameters.
%
% Inputs:
%   Wrr    - Recurrent connectivity matrix
%   Wrf    - Feedforward connectivity matrix
%   params - Network parameters
%
% Example:
%   validate_connectivity(Wrr, Wrf, params);

    % Check dimensions
    expected_rr_size = params.Ne * sum(params.Kr(:,1)) + params.Ni * sum(params.Kr(:,2));
    expected_rf_size = params.Nx * sum(params.Kx);
    
    if numel(Wrr) ~= expected_rr_size
        warning('SpikingSim:ConnectivitySize', ...
                'Recurrent connectivity size mismatch');
    end
    
    if numel(Wrf) ~= expected_rf_size
        warning('SpikingSim:ConnectivitySize', ...
                'Feedforward connectivity size mismatch');
    end
    
    % Check for valid indices
    if any(Wrr < 1) || any(Wrr > params.N)
        error('SpikingSim:InvalidIndices', ...
              'Invalid neuron indices in recurrent connectivity');
    end
    
    if any(Wrf < 1) || any(Wrf > params.N)
        error('SpikingSim:InvalidIndices', ...
              'Invalid neuron indices in feedforward connectivity');
    end
end
