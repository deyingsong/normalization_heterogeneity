function styleOptions = setPlotStyle(styleName, varargin)
% SETPLOTSTYLE Set consistent plot styling for different publication types
%
% Syntax:
%   opts = setPlotStyle(styleName)
%   opts = setPlotStyle(styleName, 'Name', Value, ...)
%
% Description:
%   Configure plot styling for different publication formats.
%
% Inputs:
%   styleName - Style preset: 'manuscript', 'presentation', 'poster', etc.
%   Name-Value pairs:
%     'Width'  - Figure width in cm
%     'Height' - Figure height in cm
%
% Outputs:
%   styleOptions - Structure with all style settings
%
% Example:
%   style = setPlotStyle('manuscript', 'Width', 8.5, 'Height', 6);
%   set(gcf, style.FigureOpt{:});

    % Constants
    DEFAULT_FONT = 'Arial';
    UNITS = 'centimeters';
    MAX_A4_HEIGHT = 23.5; % cm
    
    % Parse additional parameters
    params = utils_plot.parseNameValuePairs(varargin{:});
    params = utils_plot.validateField(params, 'Width', 8.5, 'numeric');
    params = utils_plot.validateField(params, 'Height', 6, 'numeric');
    
    % Initialize style structure
    style = struct();
    style.FontName = DEFAULT_FONT;
    style.Units = UNITS;
    
    % Define style presets
    switch lower(styleName)
        case 'manuscript'
            style.LabelSize = 12;
            style.TitleSize = 10;
            style.AxisLabelSize = 9;
            style.AxisSize = 8;
            style.LegendSize = 8;
            style.TextSize = 7;
            style.LineWidth = 1;
            style.MarkerSize = 6;
            validateHeight(params.Height, MAX_A4_HEIGHT, 'manuscript');
            
        case 'presentation'
            style.LabelSize = 14;
            style.TitleSize = 16;
            style.AxisLabelSize = 12;
            style.AxisSize = 10;
            style.LegendSize = 10;
            style.TextSize = 10;
            style.LineWidth = 2;
            style.MarkerSize = 8;
            
        case 'poster'
            style.LabelSize = 24;
            style.TitleSize = 28;
            style.AxisLabelSize = 20;
            style.AxisSize = 18;
            style.LegendSize = 18;
            style.TextSize = 16;
            style.LineWidth = 3;
            style.MarkerSize = 10;
            
        case 'thesis'
            style.LabelSize = 11;
            style.TitleSize = 12;
            style.AxisLabelSize = 10;
            style.AxisSize = 9;
            style.LegendSize = 9;
            style.TextSize = 8;
            style.LineWidth = 1;
            style.MarkerSize = 5;
            validateHeight(params.Height, MAX_A4_HEIGHT, 'thesis');
            
        case 'custom'
            % Use user-provided values or defaults
            style.LabelSize = getFieldOrDefault(params, 'LabelSize', 8);
            style.TitleSize = getFieldOrDefault(params, 'TitleSize', 9);
            style.AxisLabelSize = getFieldOrDefault(params, 'AxisLabelSize', 8);
            style.AxisSize = getFieldOrDefault(params, 'AxisSize', 8);
            style.LegendSize = getFieldOrDefault(params, 'LegendSize', 8);
            style.TextSize = getFieldOrDefault(params, 'TextSize', 7);
            style.LineWidth = getFieldOrDefault(params, 'LineWidth', 1);
            style.MarkerSize = getFieldOrDefault(params, 'MarkerSize', 6);
            
        otherwise
            error('setPlotStyle:UnknownStyle', ...
                  'Unknown style: %s. Use manuscript, presentation, poster, thesis, or custom', ...
                  styleName);
    end
    
    % Build option cell arrays
    styleOptions = buildStyleOptions(style, params);
    
    % Apply to current figure if it exists
    applyStyleToFigure(styleOptions);
    
    % --- Helper Functions ---
    function validateHeight(height, maxHeight, context)
        if height > maxHeight
            warning('setPlotStyle:HeightExceeded', ...
                    '%s figure height (%.1f cm) exceeds maximum (%.1f cm)', ...
                    context, height, maxHeight);
        end
    end
    
    function value = getFieldOrDefault(params, field, default)
        if isfield(params, field)
            value = params.(field);
        else
            value = default;
        end
    end
    
    function opts = buildStyleOptions(style, params)
        % Build option cell arrays for different elements
        opts = struct();
        
        % Figure options
        opts.FigureOpt = {
            'PaperUnits', style.Units, ...
            'PaperSize', [params.Width, params.Height], ...
            'PaperPosition', [0, 0, params.Width, params.Height], ...
            'Color', 'white'
        };
        
        % Axes options
        opts.AxisOpt = {
            'FontSize', style.AxisSize, ...
            'FontName', style.FontName, ...
            'LineWidth', style.LineWidth * 0.75, ...
            'TickDir', 'out', ...
            'Box', 'off'
        };
        
        % Label options
        opts.LabelOpt = {
            'FontSize', style.LabelSize, ...
            'FontName', style.FontName, ...
            'FontWeight', 'bold'
        };
        
        % Axis label options
        opts.AxisLabelOpt = {
            'FontSize', style.AxisLabelSize, ...
            'FontName', style.FontName
        };
        
        % Title options
        opts.TitleOpt = {
            'FontSize', style.TitleSize, ...
            'FontName', style.FontName, ...
            'FontWeight', 'bold'
        };
        
        % Legend options
        opts.LegendOpt = {
            'FontSize', style.LegendSize, ...
            'FontName', style.FontName
        };
        
        % Line options
        opts.LineOpt = {
            'LineWidth', style.LineWidth
        };
        
        % Marker options
        opts.MarkerOpt = {
            'MarkerSize', style.MarkerSize
        };
        
        % Store dimensions
        opts.Width = params.Width;
        opts.Height = params.Height;
        opts.StyleName = style;
    end
    
    function applyStyleToFigure(opts)
        try
            set(gcf, opts.FigureOpt{:});
            utils_plot.setFigureProperties('AxisOpt', opts.AxisOpt, ...
                              'LabelOpt', opts.AxisLabelOpt, ...
                              'TitleOpt', opts.TitleOpt);
        catch ME
            warning('setPlotStyle:ApplyError', ...
                    'Could not apply style: %s', ME.message);
        end
    end
end