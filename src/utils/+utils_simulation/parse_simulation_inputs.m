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
        'save', true, ...
        'CompCorr', false, ...
        'loadfr1', true, ...
        'fixW', true, ...
        'useWfile', true, ...
        'plotPopR', false, ...
        'savecurrent', false, ...
        'Layer1only', false, ...
        'saveRm', false, ...
        'saveSx', false, ...
        'saveS2', true, ...
        'saveParam', true, ...
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
