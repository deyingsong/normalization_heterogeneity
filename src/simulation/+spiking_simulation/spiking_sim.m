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
    [options, param_changes] = parse_simulation_inputs(varargin{:});
    
    % Initialize parameters
    params = initialize_network_parameters();
    
    % Apply parameter changes
    params = apply_parameter_changes(params, param_changes);
    
    % Validate parameters
    validate_parameters(params, options);
    
    % Run appropriate simulation
    switch lower(simulation_type)
        case 'default'
            run_default_simulation(params, options);
        case 'broad_weight'
            run_broad_weight_simulation(params, options);
        case 'current_noise'
            run_current_noise_simulation(params, options);
        case 'current_gaussian_noise'
            run_gaussian_noise_simulation(params, options);
        case 'fisher_info_default'
            run_fisher_info_simulation(params, options, false);
        case 'fisher_info_match_indegree'
            run_fisher_info_simulation(params, options, true);
        case 'match_indegree'
            run_match_indegree_simulation(params, options);
        otherwise
            error('SpikingSim:InvalidType', ...
                  'Unknown simulation type: %s', simulation_type);
    end
end


%% ========================================================================
%  INPUT PARSING AND VALIDATION
%% ========================================================================

function [options, param_changes] = parse_simulation_inputs(varargin)
% PARSE_SIMULATION_INPUTS - Parse and validate input arguments
%
% H1 Line: Parse optional inputs for simulation configuration
%
% Syntax:
%   [options, param_changes] = parse_simulation_inputs()
%   [options, param_changes] = parse_simulation_inputs(options)
%   [options, param_changes] = parse_simulation_inputs(options, param_changes)
%
% Description:
%   Processes variable input arguments and sets default values for
%   missing options. Ensures all required fields are present.
%
% Inputs:
%   varargin - Variable arguments: options struct and/or param_changes cell
%
% Outputs:
%   options       - Validated options structure
%   param_changes - Cell array of parameter modifications
%
% Example:
%   [options, param_changes] = parse_simulation_inputs(opt_struct, changes);

    % Initialize defaults
    options = struct();
    param_changes = {};
    
    % Parse inputs based on number of arguments
    if nargin >= 1 && isstruct(varargin{1})
        options = varargin{1};
    end
    if nargin >= 2 && iscell(varargin{2})
        param_changes = varargin{2};
    end
    
    % Set default option values
    default_options = struct(...
        'save', false, ...
        'CompCorr', false, ...
        'loadfr1', true, ...
        'fixW', false, ...
        'useWfile', false, ...
        'plotPopR', false, ...
        'savecurrent', false, ...
        'Layer1only', false, ...
        'saveRm', false, ...
        'saveSx', false, ...
        'saveS2', true, ...
        'saveParam', false, ...
        'saveW', false ...
    );
    
    % Merge with defaults
    option_fields = fieldnames(default_options);
    for i = 1:length(option_fields)
        if ~isfield(options, option_fields{i})
            options.(option_fields{i}) = default_options.(option_fields{i});
        end
    end
end

function validate_parameters(params, options)
% VALIDATE_PARAMETERS - Check parameter consistency and requirements
%
% H1 Line: Validate simulation parameters and options
%
% Syntax:
%   validate_parameters(params, options)
%
% Description:
%   Checks that all required parameters are present and consistent.
%   Throws informative errors if validation fails.
%
% Inputs:
%   params  - Parameter structure
%   options - Options structure
%
% Example:
%   validate_parameters(params, options);

    % Check required fields for saving
    if options.save && ~isfield(options, 'filename')
        error('SpikingSim:MissingFilename', ...
              'Options.save is true but no filename specified');
    end
    
    % Check correlation computation requirements
    if options.CompCorr && ~isfield(options, 'Nc')
        error('SpikingSim:MissingNc', ...
              'Options.CompCorr is true but Nc not specified');
    end
    
    % Check file loading requirements
    if options.loadfr1 && ~isfield(options, 'fr_fname1')
        error('SpikingSim:MissingFrFile', ...
              'Options.loadfr1 is true but fr_fname1 not specified');
    end
    
    if options.useWfile && ~isfield(options, 'W_fname')
        error('SpikingSim:MissingWFile', ...
              'Options.useWfile is true but W_fname not specified');
    end
    
    % Validate network size
    if params.N > 200000
        error('SpikingSim:NetworkTooLarge', ...
              'Network size N=%d exceeds maximum (200000)', params.N);
    end
    
    % Validate numerical parameters
    if params.dt <= 0
        error('SpikingSim:InvalidTimeStep', ...
              'Time step dt must be positive');
    end
    
    if params.T <= 0
        error('SpikingSim:InvalidDuration', ...
              'Simulation duration T must be positive');
    end
    
    % Check array dimensions match
    if size(params.Iapp, 2) ~= size(params.Jx, 2)
        error('SpikingSim:DimensionMismatch', ...
              'Iapp and Jx must have same number of parameter sets');
    end
end

function params = apply_parameter_changes(params, param_changes)
% APPLY_PARAMETER_CHANGES - Apply user-specified parameter modifications
%
% H1 Line: Modify default parameters with user specifications
%
% Syntax:
%   params = apply_parameter_changes(params, param_changes)
%
% Description:
%   Updates parameter structure with values specified in param_changes.
%   Validates changes and provides warnings for unknown parameters.
%
% Inputs:
%   params        - Original parameter structure
%   param_changes - Cell array with {param_name, value} pairs
%
% Outputs:
%   params - Modified parameter structure
%
% Example:
%   changes = {'T', 10000; 'dt', 0.1};
%   params = apply_parameter_changes(params, changes);

    if isempty(param_changes)
        return;
    end
    
    % Validate param_changes format
    if ~iscell(param_changes) || size(param_changes, 2) ~= 2
        error('SpikingSim:InvalidParamChanges', ...
              'param_changes must be a cell array with 2 columns');
    end
    
    % Apply each change
    for i = 1:size(param_changes, 1)
        param_name = param_changes{i, 1};
        param_value = param_changes{i, 2};
        
        % Check if parameter exists
        if ~isfield(params, param_name)
            warning('SpikingSim:UnknownParameter', ...
                    'Unknown parameter: %s', param_name);
        end
        
        % Apply the change
        params.(param_name) = param_value;
    end
    
    % Recalculate derived parameters
    params = update_derived_parameters(params);
end

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
