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
    expected_rr_size = params.N * sum(params.Kr(:));
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
