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
        
        % Apply the change
        params.(param_name) = param_value;
    end
    
    % Recalculate derived parameters
    params = utils_simulation.update_derived_parameters(params);
end
