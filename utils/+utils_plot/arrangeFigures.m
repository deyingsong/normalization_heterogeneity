function figHandles = arrangeFigures(varargin)
% ARRANGEFIGURES Arrange multiple figure windows in a grid
%
% Syntax:
%   arrangeFigures()
%   arrangeFigures('Name', Value, ...)
%   handles = arrangeFigures(...)
%
% Description:
%   Automatically arranges open figure windows in a grid layout.
%
% Inputs (Name-Value pairs):
%   'Figures'  - Figure handles or 'all' (default: 'all')
%   'Grid'     - [rows, cols] grid size or 'auto' (default: 'auto')
%   'Margin'   - Margin between figures in pixels (default: 10)
%   'Screen'   - Screen number for multi-monitor (default: 1)
%
% Outputs:
%   figHandles - Handles of arranged figures
%
% Example:
%   arrangeFigures('Grid', [2, 3], 'Margin', 20);

    % Parse parameters
    params = parseNameValuePairs(varargin{:});
    params = validateField(params, 'Figures', 'all');
    params = validateField(params, 'Grid', 'auto');
    params = validateField(params, 'Margin', 10, 'numeric');
    params = validateField(params, 'Screen', 1, 'numeric');
    
    % Get figures to arrange
    if ischar(params.Figures) && strcmpi(params.Figures, 'all')
        figHandles = findobj('Type', 'figure');
        figHandles = sort(figHandles);
    else
        figHandles = params.Figures;
    end
    
    numFigs = length(figHandles);
    if numFigs == 0
        warning('arrangeFigures:NoFigures', 'No figures to arrange');
        return;
    end
    
    % Determine grid size
    if ischar(params.Grid) && strcmpi(params.Grid, 'auto')
        gridSize = determineOptimalGrid(numFigs);
    else
        gridSize = params.Grid;
        if numel(gridSize) ~= 2
            error('arrangeFigures:InvalidGrid', 'Grid must be [rows, cols]');
        end
    end
    
    % Get screen dimensions
    screenSize = getScreenSize(params.Screen);
    
    % Calculate figure positions
    positions = calculateGridPositions(screenSize, gridSize, params.Margin);
    
    % Arrange figures
    for i = 1:min(numFigs, numel(positions))
        try
            set(figHandles(i), 'Position', positions{i});
        catch ME
            warning('arrangeFigures:PositionError', ...
                    'Could not position figure %d: %s', figHandles(i), ME.message);
        end
    end
    
    % --- Helper Functions ---
    function grid = determineOptimalGrid(n)
        % Calculate optimal grid for n figures
        cols = ceil(sqrt(n));
        rows = ceil(n / cols);
        grid = [rows, cols];
    end
    
    function screenRect = getScreenSize(screenNum)
        monitors = get(0, 'MonitorPositions');
        if screenNum > size(monitors, 1)
            warning('arrangeFigures:InvalidScreen', ...
                    'Screen %d not found, using primary screen', screenNum);
            screenNum = 1;
        end
        screenRect = monitors(screenNum, :);
        % Adjust for taskbar (approximate)
        screenRect(4) = screenRect(4) - 80;
    end
    
    function positions = calculateGridPositions(screenSize, gridSize, margin)
        % Calculate position for each grid cell
        availWidth = screenSize(3) - screenSize(1);
        availHeight = screenSize(4) - screenSize(2);
        
        figWidth = (availWidth - margin * (gridSize(2) + 1)) / gridSize(2);
        figHeight = (availHeight - margin * (gridSize(1) + 1)) / gridSize(1);
        
        positions = cell(gridSize(1), gridSize(2));
        for row = 1:gridSize(1)
            for col = 1:gridSize(2)
                x = screenSize(1) + margin + (col - 1) * (figWidth + margin);
                y = screenSize(2) + availHeight - row * (figHeight + margin);
                positions{row, col} = round([x, y, figWidth, figHeight]);
            end
        end
        positions = positions(:);
    end
end
