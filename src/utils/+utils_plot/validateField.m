function params = validateField(params, fieldName, defaultValue, dataType)
% VALIDATEFIELD Validate and set default values for structure fields
%
% Syntax:
%   params = validateField(params, fieldName, defaultValue)
%   params = validateField(params, fieldName, defaultValue, dataType)
%
% Description:
%   Validates structure fields, sets defaults, and checks data types.
%   Case-insensitive field matching with correction notifications.
%
% Inputs:
%   params       - Structure to validate
%   fieldName    - Name of field to check
%   defaultValue - Default value if field doesn't exist
%   dataType     - Expected data type (optional)
%
% Outputs:
%   params - Updated structure with validated field
%
% Example:
%   opts = struct('color', 'red');
%   opts = validateField(opts, 'LineWidth', 2, 'numeric');

    % Input validation
    if nargin < 3
        error('validateField:InsufficientArgs', ...
              'At least 3 arguments required: params, fieldName, defaultValue');
    end
    
    if ~isstruct(params) && ~isempty(params)
        error('validateField:InvalidInput', 'First argument must be a structure');
    end
    
    if ~ischar(fieldName) && ~isstring(fieldName)
        error('validateField:InvalidFieldName', 'Field name must be a string');
    end
    
    % Initialize if empty
    if isempty(params)
        params = struct();
    end
    
    % Get existing field names
    existingFields = fieldnames(params);
    
    % Check for case-insensitive match
    matchIdx = strcmpi(existingFields, fieldName);
    
    if any(matchIdx)
        % Field exists (possibly with different case)
        actualField = existingFields{matchIdx};
        
        if ~strcmp(actualField, fieldName)
            % Case mismatch - correct it
            fprintf('validateField: Correcting field name "%s" to "%s"\n', ...
                    actualField, fieldName);
            params.(fieldName) = params.(actualField);
            params = rmfield(params, actualField);
        end
        
        % Validate data type if specified
        if nargin >= 4 && ~isempty(dataType)
            validateDataType(params.(fieldName), fieldName, dataType);
        end
    else
        % Field doesn't exist - use default
        params.(fieldName) = defaultValue;
        
        % Log default usage if not empty
        if ~isempty(defaultValue)
            fprintf('validateField: Using default value for "%s"\n', fieldName);
        end
    end
end

function validateDataType(value, fieldName, expectedType)
    % Validate data type
    switch lower(expectedType)
        case 'numeric'
            if ~isnumeric(value)
                error('validateField:TypeMismatch', ...
                      'Field "%s" must be numeric', fieldName);
            end
        case 'char'
            if ~ischar(value) && ~isstring(value)
                error('validateField:TypeMismatch', ...
                      'Field "%s" must be a string', fieldName);
            end
        case 'logical'
            if ~islogical(value)
                error('validateField:TypeMismatch', ...
                      'Field "%s" must be logical', fieldName);
            end
        case 'cell'
            if ~iscell(value)
                error('validateField:TypeMismatch', ...
                      'Field "%s" must be a cell array', fieldName);
            end
    end
end
