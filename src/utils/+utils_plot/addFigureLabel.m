function handle = addFigureLabel(labelText, position, varargin)
% ADDFIGURELABEL Add a label to a figure at specified position
%
% Syntax:
%   h = addFigureLabel(labelText, position)
%   h = addFigureLabel(labelText, position, 'Name', Value, ...)
%
% Description:
%   Adds a text label to the current axes at a normalized position.
%
% Inputs:
%   labelText - Text string for the label
%   position  - [x, y] normalized position (0-1)
%   Name-Value pairs for text properties
%
% Outputs:
%   handle - Handle to the text object
%
% Example:
%   h = addFigureLabel('A', [0.05, 0.95], 'FontSize', 12, 'FontWeight', 'bold');

    % Validate inputs
    if ~ischar(labelText) && ~isstring(labelText)
        error('addFigureLabel:InvalidInput', 'Label text must be a string');
    end
    
    if ~isnumeric(position) || length(position) ~= 2
        error('addFigureLabel:InvalidPosition', 'Position must be [x, y]');
    end
    
    
    % Get current axes position for scaling
    axisPos = get(gca, 'Position');
    
    % Scale position relative to axes
    scaledPos = [position(1) / axisPos(3), ...
                 position(2) / axisPos(4) + 1];
    
    % Create text object
    handle = text(0, 0, labelText);
    
    % Set default properties
    set(handle, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'left', ...
        'Position', scaledPos, ...
        'Tag', 'FigureLabel', ...
        varargin{:});
end
