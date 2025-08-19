%% Main Simulation Functions Suite
% Refactored spiking neural network simulation code
% Author: Deying Song
% Date: Aug 15th, 2025
% Version: 2.0

%% ========================================================================
%  MAIN SIMULATION WRAPPER
%% ========================================================================

function spiking_sim(simulation_type, varargin)
% SPIKING_SIM - Main wrapper for spiking neural network simulations
%
% H1 Line: Runs spiking neural network simulations with various configurations
%
% Syntax:
%   spiking_sim('default', options, param_changes)
%   spiking_sim('broad_weight', options, param_changes)
%   spiking_sim('current_noise', options, param_changes)
%   spiking_sim('fisher_info', options, param_changes)
%
% Description:
%   Main entry point for running spiking neural network simulations.
%   Supports multiple simulation types with different noise models and
%   connectivity patterns.
%
% Inputs:
%   simulation_type - String specifying simulation variant:
%                    'default', 'broad_weight', 'current_noise',
%                    'current_gaussian_noise', 'fisher_info_default',
%                    'fisher_info_match_indegree', 'match_indegree'
%   options        - Structure with simulation options (optional)
%   param_changes  - Cell array of parameter modifications (optional)
%
% Outputs:
%   Results saved to specified filename in options
%
% Example:
%   options = struct('save', 1, 'filename', 'results.mat');
%   param_changes = {'T', 10000; 'dt', 0.05};
%   spiking_sim('default', options, param_changes);

    % Parse inputs
    [options, param_changes] = utils_simulation.parse_simulation_inputs(varargin{:});
    
    % Initialize parameters
    params = utils_simulation.initialize_network_parameters();
    
    % Apply parameter changes
    params = utils_simulation.apply_parameter_changes(params, param_changes);
    
    % Validate parameters
    utils_simulation.validate_parameters(params, options);
    
    % Run appropriate simulation
    switch lower(simulation_type)
        case 'default'
            spiking_simulation.run_default_simulation(params, options);
        case 'broad_weight'
            spiking_simulation.run_broad_weight_simulation(params, options);
        case 'current_noise'
            spiking_simulation.run_current_noise_simulation(params, options);
        case 'current_gaussian_noise'
            spiking_simulation.run_gaussian_noise_simulation(params, options);
        case 'fisher_info_default'
            spiking_simulation.run_fisher_info_simulation(params, options, false);
        case 'fisher_info_match_indegree'
            spiking_simulation.run_fisher_info_simulation(params, options, true);
        case 'match_indegree'
            spiking_simulation.run_match_indegree_simulation(params, options);
        otherwise
            error('SpikingSim:InvalidType', ...
                  'Unknown simulation type: %s', simulation_type);
    end
end


