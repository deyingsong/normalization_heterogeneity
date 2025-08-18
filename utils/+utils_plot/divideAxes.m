function [positions, axesHandles] = divideAxes(varargin)
% DIVIDEAXES Divide figure into custom grid of axes with specified spacing
%
% Syntax:
%   positions = divideAxes(divX, divY)
%   positions = divideAxes(divX, divY, bounds)
%   positions = divideAxes(divX, divY, bounds, sepX, sepY)
%   [positions, axesHandles] = divideAxes(..., 'create')
%
% Description:
%   Creates a grid of axes positions with custom relative sizes and spacing.
%   Can optionally create the axes automatically.
%
% Inputs:
%   divX    - Vector of relative sizes in X direction (default: [1,1])
%   divY    - Vector of relative sizes in Y direction (default: [1,1])
%   bounds  - [x0, y0, width, height] outer boundaries (default: [0.1,0.1,0.75,0.75])
%   sepX    - Separations between columns (length = length(divX)-1)
%   sepY    - Separations between rows (length = length(divY)-1)
%   option  - 'create' to automatically create axes
%
% Outputs:
%   positions    - Cell array of axis positions [x,y,width,height]
%   axesHandles  - Array of axes handles (if 'create' option used)
%
% Example:
%   % Create 2x3 grid with different column widths
%   [pos, ax] = divideAxes([1,2,1], [1,1], [0.1,0.1,0.8,0.8], 0.05, 0.05, 'create');

    % Constants
    DEFAULT_BOUNDS = [0.1, 0.85, 0.1, 0.85];  % [x0, x_width, y0, y_height]
    DEFAULT_SEP_X = 0.04;  % Default X separation (normalized)
    DEFAULT_SEP_Y = 0.07;  % Default Y separation (normalized)
    
    % Parse inputs
    [divX, divY, xBounds, yBounds, sepX, sepY, createAxes] = parseInputs(varargin{:});
    
    % Validate inputs
    validateDivisions(divX, divY);
    validateSeparations(sepX, sepY, divX, divY);
    
    % Handle scalar divisions (convert to equal divisions)
    divX = expandScalarDivision(divX);
    divY = expandScalarDivision(divY);
    
    % Expand scalar separations to vectors
    sepX = expandScalarSeparation(sepX, length(divX) - 1);
    sepY = expandScalarSeparation(sepY, length(divY) - 1);
    
    % Calculate normalized sizes
    [normDivX, normSepX] = normalizeValues(divX, sepX, xBounds(2));
    [normDivY, normSepY] = normalizeValues(divY, sepY, yBounds(2));
    
    % Flip Y for MATLAB coordinate system
    normDivY = fliplr(normDivY);
    normSepY = fliplr(normSepY);
    
    % Calculate positions
    positions = calculatePositions(normDivX, normDivY, normSepX, normSepY, xBounds(1), yBounds(1));
    positions = flipud(positions);
    
    % Create axes if requested
    axesHandles = [];
    if createAxes
        axesHandles = createAxesFromPositions(positions);
    end
    
    % --- Nested Functions ---
    function [dx, dy, xb, yb, sx, sy, create] = parseInputs(varargin)
        % Set defaults
        dx = [1, 1];
        dy = [1, 1];
        xb = DEFAULT_BOUNDS([1,2]);
        yb = DEFAULT_BOUNDS([3,4]);
        sx = DEFAULT_SEP_X;
        sy = DEFAULT_SEP_Y;
        create = false;
        
        % Check for 'create' option
        if nargin > 0 && ischar(varargin{end})
            if strcmpi(varargin{end}, 'create')
                create = true;
                varargin = varargin(1:end-1);
            end
        end
        
        % Parse positional arguments
        if length(varargin) >= 1 && ~isempty(varargin{1})
            dx = varargin{1};
        end
        if length(varargin) >= 2 && ~isempty(varargin{2})
            dy = varargin{2};
        end
        if length(varargin) >= 3 && ~isempty(varargin{3})
            bounds = varargin{3};
            if length(bounds) == 4
                xb = bounds([1,3]);
                yb = bounds([2,4]);
            else
                error('divideAxes:InvalidBounds', 'Bounds must be [x0, y0, width, height]');
            end
        end
        if length(varargin) >= 4 && ~isempty(varargin{4})
            sx = varargin{4};
        end
        if length(varargin) >= 5 && ~isempty(varargin{5})
            sy = varargin{5};
        end
    end
    
    function validateDivisions(dx, dy)
        if ~isnumeric(dx) || ~isnumeric(dy)
            error('divideAxes:InvalidInput', 'Divisions must be numeric');
        end
        if any(dx <= 0) || any(dy <= 0)
            error('divideAxes:InvalidDivision', 'All divisions must be positive');
        end
    end
    
    function validateSeparations(sx, sy, dx, dy)
        if ~isnumeric(sx) || ~isnumeric(sy)
            error('divideAxes:InvalidSeparation', 'Separations must be numeric');
        end
        if any(sx < 0) || any(sy < 0)
            error('divideAxes:NegativeSeparation', 'Separations cannot be negative');
        end
    end
    
    function div = expandScalarDivision(div)
        if length(div) == 1 && div > 1 && div == round(div)
            div = ones(1, div);
        end
    end
    
    function sep = expandScalarSeparation(sep, targetLength)
        if length(sep) == 1 && targetLength > 0
            sep = repmat(sep, 1, targetLength);
        end
    end
    
    function [normDiv, normSep] = normalizeValues(div, sep, scale)
        total = sum(div) + sum(sep);
        normDiv = (div / total) * scale;
        normSep = (sep / total) * scale;
    end
    
    function pos = calculatePositions(dx, dy, sx, sy, x0, y0)
        nRows = length(dy);
        nCols = length(dx);
        pos = cell(nRows, nCols);
        
        for i = 1:nRows
            for j = 1:nCols
                xPos = x0 + sum(dx(1:j-1)) + sum(sx(1:j-1));
                yPos = y0 + sum(dy(1:i-1)) + sum(sy(1:i-1));
                pos{i,j} = [xPos, yPos, dx(j), dy(i)];
            end
        end
    end
    
    function ax = createAxesFromPositions(pos)
        ax = zeros(size(pos));
        for i = 1:numel(pos)
            ax(i) = axes('Position', pos{i});
        end
    end
end
