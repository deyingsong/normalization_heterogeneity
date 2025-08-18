function setFigureProperties(varargin)
% SETFIGUREPROPERTIES Set consistent properties across figure elements
%
% Syntax:
%   setFigureProperties()
%   setFigureProperties('Name', Value, ...)
%
% Description:
%   Apply consistent formatting to all axes, labels, and text in figure.
%
% Inputs (Name-Value pairs):
%   'Figure'     - Figure handle (default: gcf)
%   'AxisOpt'    - Cell array of axes properties
%   'LabelOpt'   - Cell array of label properties
%   'TitleOpt'   - Cell array of title properties
%   'TickSize'   - Tick length (default: 0.01)
%
% Example:
%   setFigureProperties('AxisOpt', {'FontSize', 10}, 'LabelOpt', {'FontSize', 12});

    % Parse parameters
    params = utils_plot.parseNameValuePairs(varargin{:});
    params = utils_plot.validateField(params, 'Figure', gcf);
    params = utils_plot.validateField(params, 'AxisOpt', {}, 'cell');
    params = utils_plot.validateField(params, 'LabelOpt', {}, 'cell');
    params = utils_plot.validateField(params, 'TitleOpt', {}, 'cell');
    params = utils_plot.validateField(params, 'TickSize', 0.01, 'numeric');
    
    % Get all axes in figure
    allAxes = findobj(params.Figure, 'Type', 'axes');
    
    % Apply axes properties
    if ~isempty(params.AxisOpt)
        try
            set(allAxes, params.AxisOpt{:});
        catch ME
            warning('setFigureProperties:AxisError', ...
                    'Could not set axes properties: %s', ME.message);
        end
    end
    
    % Apply to each axes individually
    for ax = allAxes'
        % Apply label properties
        if ~isempty(params.LabelOpt)
            labels = [get(ax, 'XLabel'), get(ax, 'YLabel'), get(ax, 'ZLabel')];
            try
                set(labels, params.LabelOpt{:});
            catch ME
                warning('setFigureProperties:LabelError', ...
                        'Could not set label properties: %s', ME.message);
            end
        end
        
        % Apply title properties
        if ~isempty(params.TitleOpt)
            titleHandle = get(ax, 'Title');
            try
                set(titleHandle, params.TitleOpt{:});
            catch ME
                warning('setFigureProperties:TitleError', ...
                        'Could not set title properties: %s', ME.message);
            end
        end
        
        % Set tick size
        set(ax, 'TickLength', [params.TickSize, params.TickSize]);
    end
end
