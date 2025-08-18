function handles = plotErrorArea(x, mean_y, error_y, varargin)
% PLOTERRORAREA Plot line with shaded error region
%
% Syntax:
%   h = plotErrorArea(x, mean_y, error_y)
%   h = plotErrorArea(x, mean_y, error_y, 'Name', Value, ...)
%
% Description:
%   Creates a plot with a mean line and shaded error regions.
%   Handles NaN values and asymmetric error bounds.
%
% Inputs:
%   x       - X values (vector)
%   mean_y  - Mean Y values (vector)
%   error_y - Error values (vector or 2-column matrix for asymmetric)
%   Name-Value pairs:
%     'Color'      - Line color (default: [0,0,1])
%     'FaceColor'  - Error region color (default: lighter line color)
%     'Alpha'      - Transparency (0-1, default: 0.5)
%     'LineWidth'  - Line width (default: 1)
%     'LineStyle'  - Line style (default: '-')
%
% Outputs:
%   handles - [lineHandle, patchHandle]
%
% Example:
%   x = 1:100;
%   y = sin(x/10);
%   err = 0.1 * rand(size(x));
%   h = plotErrorArea(x, y, err, 'Color', 'r', 'Alpha', 0.3);

    % Parse input arguments
    params = parseNameValuePairs(varargin{:});
    
    % Set defaults
    params = validateField(params, 'Color', [0, 0, 1], 'numeric');
    params = validateField(params, 'FaceColor', [], 'numeric');
    params = validateField(params, 'Alpha', 0.5, 'numeric');
    params = validateField(params, 'LineWidth', 1, 'numeric');
    params = validateField(params, 'LineStyle', '-', 'char');
    
    % Auto-generate face color if not specified
    if isempty(params.FaceColor)
        params.FaceColor = lightenColor(params.Color, 0.5);
    end
    
    % Validate inputs
    validateErrorAreaInputs(x, mean_y, error_y);
    
    % Ensure column vectors
    x = x(:);
    mean_y = mean_y(:);
    
    % Process error values
    error_y = processErrorValues(error_y, length(mean_y));
    
    % Remove NaN values
    [x, mean_y, error_y] = removeNaNValues(x, mean_y, error_y);
    
    % Check if data is valid after NaN removal
    if isempty(x)
        warning('plotErrorArea:NoValidData', 'No valid data points after removing NaNs');
        handles = [NaN, NaN];
        return;
    end
    
    % Set renderer for transparency
    if params.Alpha < 1
        set(gcf, 'Renderer', 'OpenGL');
    end
    
    % Create error patch
    patchHandle = createErrorPatch(x, mean_y, error_y, params);
    hold on;
    
    % Create mean line
    lineHandle = [];
    if ~strcmp(params.LineStyle, 'none')
        lineHandle = plot(x, mean_y, ...
                         'Color', params.Color, ...
                         'LineWidth', params.LineWidth, ...
                         'LineStyle', params.LineStyle);
    end
    
    % Return handles
    handles = [lineHandle, patchHandle];
    
    % --- Helper Functions ---
    function validateErrorAreaInputs(x, y, err)
        if ~isvector(x) || ~isvector(y)
            error('plotErrorArea:InvalidInput', 'X and Y must be vectors');
        end
        if length(x) ~= length(y)
            error('plotErrorArea:SizeMismatch', 'X and Y must have the same length');
        end
        if size(err, 1) ~= length(y) && size(err, 2) ~= length(y)
            error('plotErrorArea:ErrorSizeMismatch', ...
                  'Error must match the length of data');
        end
    end
    
    function err = processErrorValues(err, dataLength)
        % Ensure error is properly oriented
        if size(err, 1) ~= dataLength && size(err, 2) == dataLength
            err = err';
        end
        
        % Handle single column (symmetric errors)
        if size(err, 2) == 1
            err(:, 2) = err(:, 1);
        end
    end
    
    function [xClean, yClean, errClean] = removeNaNValues(x, y, err)
        nanIdx = isnan(y) | any(isnan(err), 2);
        validIdx = ~nanIdx;
        xClean = x(validIdx);
        yClean = y(validIdx);
        errClean = err(validIdx, :);
    end
    
    function h = createErrorPatch(x, y, err, params)
        % Create vertices for patch
        xPatch = [x; flipud(x)];
        yPatch = [y - err(:,1); flipud(y + err(:,2))];
        
        h = patch(xPatch, yPatch, params.FaceColor, ...
                 'FaceAlpha', params.Alpha, ...
                 'EdgeColor', 'none');
    end
    
    function lightColor = lightenColor(color, factor)
        % Lighten a color by mixing with white
        white = [1, 1, 1];
        lightColor = color + (white - color) * factor;
        lightColor = min(lightColor, 1); % Clamp to [0,1]
    end
end
