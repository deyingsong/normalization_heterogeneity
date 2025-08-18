function handles = plotFilledArea(x, y, varargin)
% PLOTFILLEDAREA Plot filled area under curve
%
% Syntax:
%   h = plotFilledArea(x, y)
%   h = plotFilledArea(x, y, 'Name', Value, ...)
%
% Description:
%   Creates a filled area plot from baseline to curve.
%
% Inputs:
%   x - X values (vector)
%   y - Y values (vector)
%   Name-Value pairs:
%     'Color'      - Line color (default: [0,0,1])
%     'FaceColor'  - Fill color (default: lighter line color)
%     'Alpha'      - Transparency (default: 0.5)
%     'Baseline'   - Y value for baseline (default: 0)
%
% Outputs:
%   handles - [lineHandle, patchHandle]
%
% Example:
%   x = linspace(0, 2*pi, 100);
%   h = plotFilledArea(x, sin(x), 'Color', 'b', 'Alpha', 0.3);

    % Parse parameters
    params = parseNameValuePairs(varargin{:});
    params = validateField(params, 'Color', [0, 0, 1], 'numeric');
    params = validateField(params, 'FaceColor', [], 'numeric');
    params = validateField(params, 'Alpha', 0.5, 'numeric');
    params = validateField(params, 'LineWidth', 1, 'numeric');
    params = validateField(params, 'LineStyle', '-', 'char');
    params = validateField(params, 'Baseline', 0, 'numeric');
    
    % Auto-generate face color
    if isempty(params.FaceColor)
        params.FaceColor = lightenColor(params.Color, 0.5);
    end
    
    % Validate inputs
    if ~isvector(x) || ~isvector(y)
        error('plotFilledArea:InvalidInput', 'X and Y must be vectors');
    end
    if length(x) ~= length(y)
        error('plotFilledArea:SizeMismatch', 'X and Y must have same length');
    end
    
    % Ensure column vectors
    x = x(:);
    y = y(:);
    
    % Remove NaN values
    validIdx = ~(isnan(x) | isnan(y));
    x = x(validIdx);
    y = y(validIdx);
    
    if isempty(x)
        warning('plotFilledArea:NoData', 'No valid data after removing NaNs');
        handles = [NaN, NaN];
        return;
    end
    
    % Set renderer for transparency
    if params.Alpha < 1
        set(gcf, 'Renderer', 'OpenGL');
    end
    
    % Create filled area
    baseline = repmat(params.Baseline, size(y));
    xPatch = [x; flipud(x)];
    yPatch = [baseline; flipud(y)];
    
    patchHandle = patch(xPatch, yPatch, params.FaceColor, ...
                       'FaceAlpha', params.Alpha, ...
                       'EdgeColor', 'none');
    hold on;
    
    % Create outline
    lineHandle = [];
    if ~strcmp(params.LineStyle, 'none')
        lineHandle = plot(x, y, ...
                         'Color', params.Color, ...
                         'LineWidth', params.LineWidth, ...
                         'LineStyle', params.LineStyle);
    end
    
    handles = [lineHandle, patchHandle];
    
    % --- Helper Function ---
    function lightColor = lightenColor(color, factor)
        white = [1, 1, 1];
        lightColor = color + (white - color) * factor;
        lightColor = min(lightColor, 1);
    end
end
