function freezeColors(varargin)
%FREEZECOLORS Freeze current colormap to prevent changes
%
% Purpose:
%   Locks the current colormap so that subsequent colormap changes don't
%   affect the current figure or axes.
%
% Usage:
%   freezeColors()          % Freeze current axes
%   freezeColors(h)         % Freeze specific axes handle h

if nargin == 0
    h = gca;
else
    h = varargin{1};
end

% Get current colormap and color data
cmap = colormap(h);
clims = clim(h);

% Convert indexed color data to RGB
children = get(h, 'Children');
for i = 1:length(children)
    child = children(i);
    if isprop(child, 'CData') && ~isempty(get(child, 'CData'))
        cdata = get(child, 'CData');
        if size(cdata, 3) == 1 % Indexed color data
            % Normalize to colormap indices
            cdata_norm = (cdata - clims(1)) / (clims(2) - clims(1));
            cdata_norm = max(0, min(1, cdata_norm));
            
            % Convert to RGB
            indices = round(cdata_norm * (size(cmap, 1) - 1)) + 1;
            rgb_data = reshape(cmap(indices, :), [size(cdata), 3]);
            
            % Set RGB data
            set(child, 'CData', rgb_data);
        end
    end
end
end