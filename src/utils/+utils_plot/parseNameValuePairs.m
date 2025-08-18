function params = parseNameValuePairs(varargin)
% PARSENAMEVALUEPAIRS Parse name-value pair arguments into structure
%
% Syntax:
%   params = parseNameValuePairs('name1', value1, 'name2', value2, ...)
%
% Description:
%   Converts name-value pairs into a structure for easy access.
%
% Inputs:
%   varargin - Name-value pairs
%
% Outputs:
%   params - Structure with fields corresponding to names
%
% Example:
%   opts = parseNameValuePairs('Color', 'red', 'LineWidth', 2);

    params = struct();
    
    % Check for even number of arguments
    if mod(length(varargin), 2) ~= 0
        error('parseNameValuePairs:OddNumber', ...
              'Arguments must come in name-value pairs');
    end
    
    % Process each pair
    for i = 1:2:length(varargin)
        name = varargin{i};
        value = varargin{i+1};
        
        % Validate name
        if ~ischar(name) && ~isstring(name)
            error('parseNameValuePairs:InvalidName', ...
                  'Parameter names must be strings (argument %d)', i);
        end
        
        % Convert to valid field name
        fieldName = matlab.lang.makeValidName(char(name));
        
        % Store in structure
        params.(fieldName) = value;
    end
end
