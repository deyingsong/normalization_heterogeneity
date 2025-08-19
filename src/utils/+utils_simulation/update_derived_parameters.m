function params = update_derived_parameters(params)
% UPDATE_DERIVED_PARAMETERS - Recalculate parameters that depend on others
%
% H1 Line: Update parameters derived from base parameters
%
% Syntax:
%   params = update_derived_parameters(params)
%
% Description:
%   Recalculates parameters that depend on other parameters after
%   changes have been applied. Ensures consistency.
%
% Inputs:
%   params - Parameter structure with base parameters
%
% Outputs:
%   params - Updated parameter structure
%
% Example:
%   params = update_derived_parameters(params);

    % Update population sizes if dimensions changed
    params.Ne = params.Ne1 * params.Ne1 / 2;
    params.Ni = params.Ni1 * params.Ni1 / 2;
    params.Nx = params.Nx1 * params.Nx1 / 2;
    params.N = params.Ne + params.Ni;
    
    % Update maximum spike count
    params.maxns = params.N * params.T * params.maxrate;
    
    % Scale connection strengths by sqrt(N)
    params.Jr_scaled = params.Jr / sqrt(params.N);
    params.Jx_scaled = params.Jx / sqrt(params.N);
    
    % Calculate out-degrees
    params.Kr = ceil(params.Prr .* [params.Ne, params.Ne; 
                                    params.Ni, params.Ni]);
    params.Kx = ceil(params.Prx .* [params.Ne; params.Ni]);
    
    % Update leak voltage based on applied current
    params.vl = params.Iapp' .* [15, 10] - 60;
    
    % Set voltage initialization range
    params.V0min = params.vre(1);
    params.V0max = params.vT(1);
end
