function handle = createNarrowColorbar(location)
% CREATENARROWCOLORBAR Create a narrow colorbar
%
% Syntax:
%   h = createNarrowColorbar()
%   h = createNarrowColorbar(location)
%
% Description:
%   Creates a narrow colorbar that takes up less space than default.
%
% Inputs:
%   location - 'vert' (default) or 'horiz'
%
% Outputs:
%   handle - Handle to colorbar axes
%
% Example:
%   imagesc(peaks);
%   h = createNarrowColorbar('vert');

    % Constants
    STRIPE_WIDTH = 0.025;  % Relative width of colorbar
    EDGE_MARGIN = 0.02;    % Edge margin
    SPACE_3D = 0.1;        % Extra space for 3D plots
    SPACE_2D = 0.05;       % Space for 2D plots
    
    % Default location
    if nargin < 1
        location = 'vert';
    end
    
    % Validate input
    if ~ischar(location)
        error('createNarrowColorbar:InvalidInput', 'Location must be a string');
    end
    
    % Get current axes
    currentAxes = gca;
    
    % Determine color limits
    colorLimits = determineColorLimits(currentAxes);
    
    % Check for existing colorbar
    existingBar = findExistingColorbar(currentAxes);
    
    if ~isempty(existingBar)
        % Update existing colorbar
        axes(existingBar);
        handle = existingBar;
    else
        % Create new colorbar
        handle = createNewColorbar(currentAxes, location, colorLimits, ...
                                  STRIPE_WIDTH, EDGE_MARGIN, SPACE_2D, SPACE_3D);
    end
    
    % Return focus to original axes
    axes(currentAxes);
    
    % --- Helper Functions ---
    function limits = determineColorLimits(ax)
        % Check for images or surfaces
        children = get(ax, 'Children');
        hasImage = false;
        
        for i = 1:length(children)
            type = get(children(i), 'Type');
            if strcmp(type, 'image') || strcmp(type, 'surface')
                hasImage = true;
                mapping = get(children(i), 'CDataMapping');
                if strcmp(mapping, 'scaled')
                    limits = caxis;
                    return;
                end
            end
        end
        
        if hasImage
            limits = [1, size(colormap, 1)];
        else
            limits = caxis;
        end
    end
    
    function cb = findExistingColorbar(ax)
        % Search for existing colorbar associated with axes
        cb = [];
        images = findobj(gcf, 'Type', 'image', 'Tag', 'Colorbar');
        for i = 1:length(images)
            parent = get(images(i), 'Parent');
            userData = get(parent, 'UserData');
            if isfield(userData, 'PlotHandle') && userData.PlotHandle == ax
                cb = parent;
                return;
            end
        end
    end
    
    function cb = createNewColorbar(ax, loc, limits, stripe, edge, space2d, space3d)
        % Store original position
        originalUnits = get(ax, 'Units');
        set(ax, 'Units', 'normalized');
        originalPos = get(ax, 'Position');
        
        % Determine space based on view
        [azimuth, elevation] = view(ax);
        if azimuth == 0 && elevation == 90
            space = space2d;
        else
            space = space3d;
        end
        
        % Create colorbar based on location
        if strncmpi(loc, 'v', 1)
            % Vertical colorbar
            newAxPos = [originalPos(1), originalPos(2), ...
                       originalPos(3) * (1 - stripe - edge - space), originalPos(4)];
            cbPos = [originalPos(1) + (1 - stripe - edge) * originalPos(3), ...
                    originalPos(2), stripe * originalPos(3), originalPos(4)];
            
            set(ax, 'Position', newAxPos);
            cb = axes('Position', cbPos);
            
            % Create color image
            n = size(colormap, 1);
            image([0 1], limits, (1:n)', 'Tag', 'Colorbar');
            set(cb, 'YDir', 'normal', 'XTick', [], 'YAxisLocation', 'right');
            
        else
            % Horizontal colorbar
            newAxPos = [originalPos(1), ...
                       originalPos(2) + (stripe + space) * originalPos(4), ...
                       originalPos(3), ...
                       (1 - stripe - space) * originalPos(4)];
            cbPos = [originalPos(1), originalPos(2), ...
                    originalPos(3), stripe * originalPos(4)];
            
            set(ax, 'Position', newAxPos);
            cb = axes('Position', cbPos);
            
            % Create color image
            n = size(colormap, 1);
            image(limits, [0 1], (1:n), 'Tag', 'Colorbar');
            set(cb, 'YDir', 'normal', 'YTick', []);
        end
        
        % Store user data
        userData.PlotHandle = ax;
        userData.OriginalPosition = originalPos;
        set(cb, 'UserData', userData, 'Tag', 'ColorbarAxes');
        
        % Restore units
        set(ax, 'Units', originalUnits);
    end
end
